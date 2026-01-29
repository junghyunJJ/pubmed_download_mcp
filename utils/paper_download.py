"""
Paper download functionality using Jina Reader API and NCBI efetch
with improved rate limiting and error handling
"""

import os
import asyncio
import random
from typing import Dict, Any, Optional
import requests
from Bio import Entrez
from bs4 import BeautifulSoup


def parse_jats_xml_to_text(xml_content: bytes) -> str:
    """
    Parse JATS XML (PMC full text) to plain text for RAG.

    Extracts title, abstract, and body sections from PMC's JATS XML format
    and formats them as clean, readable text.

    Args:
        xml_content: Raw XML bytes from NCBI efetch

    Returns:
        Formatted plain text suitable for RAG and analysis
    """
    soup = BeautifulSoup(xml_content, "xml")

    sections = []

    # Extract article title
    title_tag = soup.find("article-title")
    if title_tag:
        title_text = title_tag.get_text(strip=True)
        sections.append(f"Title: {title_text}")

    # Extract authors
    authors = []
    for contrib in soup.find_all("contrib", {"contrib-type": "author"}):
        surname = contrib.find("surname")
        given_names = contrib.find("given-names")
        if surname:
            name = surname.get_text(strip=True)
            if given_names:
                name = f"{given_names.get_text(strip=True)} {name}"
            authors.append(name)
    if authors:
        sections.append(f"Authors: {', '.join(authors)}")

    # Extract journal info
    journal_title = soup.find("journal-title")
    if journal_title:
        sections.append(f"Journal: {journal_title.get_text(strip=True)}")

    # Extract publication year
    pub_date = soup.find("pub-date")
    if pub_date:
        year = pub_date.find("year")
        if year:
            sections.append(f"Year: {year.get_text(strip=True)}")

    # Extract PMC ID
    article_id = soup.find("article-id", {"pub-id-type": "pmc"})
    if article_id:
        sections.append(f"PMC ID: PMC{article_id.get_text(strip=True)}")

    # Extract DOI
    doi_tag = soup.find("article-id", {"pub-id-type": "doi"})
    if doi_tag:
        sections.append(f"DOI: {doi_tag.get_text(strip=True)}")

    sections.append("")  # Empty line before abstract

    # Extract abstract
    abstract = soup.find("abstract")
    if abstract:
        sections.append("Abstract:")
        # Handle structured abstracts with sections
        abstract_sections = abstract.find_all("sec")
        if abstract_sections:
            for sec in abstract_sections:
                sec_title = sec.find("title")
                if sec_title:
                    sections.append(f"\n{sec_title.get_text(strip=True)}:")
                for p in sec.find_all("p"):
                    sections.append(p.get_text(strip=True))
        else:
            # Simple abstract
            for p in abstract.find_all("p"):
                sections.append(p.get_text(strip=True))
        sections.append("")

    # Extract body content
    body = soup.find("body")
    if body:
        sections.append("Full Text:")
        sections.append("")

        # Process each section
        for sec in body.find_all("sec", recursive=False):
            _process_section(sec, sections, level=1)

        # If no sections, just get paragraphs
        if not body.find_all("sec", recursive=False):
            for p in body.find_all("p", recursive=False):
                sections.append(p.get_text(strip=True))
                sections.append("")

    # References section excluded for RAG efficiency
    # (reference text adds noise without adding semantic value)

    return "\n".join(sections)


def _process_section(sec, sections: list, level: int = 1):
    """Helper function to recursively process JATS sections."""
    # Get section title
    title = sec.find("title", recursive=False)
    if title:
        prefix = "#" * min(level + 1, 4)  # Markdown-style headers
        sections.append(f"{prefix} {title.get_text(strip=True)}")
        sections.append("")

    # Process paragraphs in this section
    for p in sec.find_all("p", recursive=False):
        text = p.get_text(strip=True)
        if text:
            sections.append(text)
            sections.append("")

    # Process nested sections
    for nested_sec in sec.find_all("sec", recursive=False):
        _process_section(nested_sec, sections, level + 1)


async def fetch_pmc_fulltext(
    pmcid: str,
    email: str,
    api_key: Optional[str] = None
) -> Dict[str, Any]:
    """
    Fetch full text from PMC using NCBI efetch API.

    This function provides a direct download from NCBI's PMC database,
    bypassing Jina Reader which may be blocked by PMC's DDoS protection.

    Args:
        pmcid: PMC ID (e.g., "PMC1234567" or "1234567")
        email: Email for NCBI Entrez API (required by NCBI policy)
        api_key: Optional NCBI API key for higher rate limits

    Returns:
        Dictionary with:
        - success: bool indicating if fetch was successful
        - content: str containing formatted full text
        - message: str with status message
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # Normalize PMC ID (remove "PMC" prefix if present)
    pmc_number = pmcid.replace("PMC", "").strip()

    try:
        def _fetch():
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_number,
                rettype="full",
                retmode="xml"
            )
            xml_content = handle.read()
            handle.close()
            return xml_content

        # Rate limiting for NCBI compliance
        await asyncio.sleep(0.5)
        xml_content = await asyncio.to_thread(_fetch)

        # Parse XML to text
        text_content = parse_jats_xml_to_text(xml_content)

        if not text_content or len(text_content.strip()) < 100:
            return {
                "success": False,
                "content": "",
                "message": f"PMC efetch returned empty or minimal content for {pmcid}"
            }

        return {
            "success": True,
            "content": text_content,
            "message": f"Successfully fetched full text from PMC via efetch for {pmcid}"
        }

    except Exception as e:
        return {
            "success": False,
            "content": "",
            "message": f"Error fetching PMC full text via efetch: {str(e)}"
        }


async def check_pmc_availability(pmid: str, email: str, api_key: Optional[str] = None) -> Optional[str]:
    """
    Check if a paper is available in PubMed Central using NCBI ID Converter API.

    Args:
        pmid: PubMed ID
        email: Email for Entrez API (for User-Agent compliance)
        api_key: Optional NCBI API key for higher rate limits

    Returns:
        PMC ID (e.g., "PMC1234567") if available, None otherwise
    """
    try:
        # Use NCBI ID Converter API for accurate PMC ID lookup
        def _check():
            url = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
            params = {
                "ids": pmid,
                "format": "json"
            }

            # Add API key if provided
            if api_key:
                params["api_key"] = api_key

            headers = {
                "User-Agent": f"PubMedDownloadMCP/1.0 ({email})"
            }

            response = requests.get(url, params=params, headers=headers, timeout=10)

            if response.status_code == 200:
                data = response.json()
                # Check if records exist and have pmcid
                if "records" in data and data["records"]:
                    record = data["records"][0]
                    # Only return if pmcid exists and status is not "error"
                    if "pmcid" in record and record.get("status") != "error":
                        return record["pmcid"]

            return None

        return await asyncio.to_thread(_check)
    except Exception as e:
        print(f"Error checking PMC availability: {str(e)}")
        return None


def check_cached_paper(doi: str, output_path: Optional[str] = None) -> Optional[str]:
    """
    Check if paper is already downloaded.

    Args:
        doi: Paper DOI
        output_path: Optional custom output path

    Returns:
        File path if paper exists, None otherwise
    """
    if output_path and os.path.exists(output_path):
        return output_path

    # Check default location
    sanitized_doi = doi.replace("/", "_").replace(":", "_")
    default_path = f"./papers/{sanitized_doi}.txt"
    if os.path.exists(default_path):
        return default_path

    return None


async def fetch_pubmed_abstract(pmid: str, email: str) -> Dict[str, Any]:
    """
    Fetch abstract and metadata from PubMed for papers without PMC or DOI.

    This function provides a fallback for papers that are not available in PMC
    and do not have a DOI. It fetches the abstract and metadata using the
    Entrez API and formats them as clean text.

    Args:
        pmid: PubMed ID
        email: Email for NCBI Entrez API (required by NCBI policy)

    Returns:
        Dictionary with:
        - success: bool indicating if fetch was successful
        - content: str containing formatted abstract and metadata
        - message: str with status message
    """
    Entrez.email = email

    try:
        # Fetch paper details using Entrez.efetch
        def _fetch():
            handle = Entrez.efetch(
                db="pubmed",
                id=pmid,
                retmode="xml"
            )
            records = Entrez.read(handle)
            handle.close()
            return records

        # Rate limiting for NCBI compliance
        await asyncio.sleep(0.5)
        records = await asyncio.to_thread(_fetch)

        pubmed_articles = records.get("PubmedArticle", [])
        if not pubmed_articles:
            return {
                "success": False,
                "content": "",
                "message": f"No article found for PMID {pmid}"
            }

        # Extract first article (should only be one for a single PMID)
        record = pubmed_articles[0]
        medline_citation = record.get("MedlineCitation", {})
        article = medline_citation.get("Article", {})

        # Extract title
        title = article.get("ArticleTitle", "Unknown Title")

        # Extract abstract
        abstract = ""
        if "Abstract" in article:
            abstract_texts = article["Abstract"].get("AbstractText", [])
            if isinstance(abstract_texts, list):
                abstract = " ".join(str(text) for text in abstract_texts)
            else:
                abstract = str(abstract_texts)

        if not abstract:
            abstract = "[No abstract available]"

        # Extract authors
        authors = []
        if "AuthorList" in article:
            for author in article["AuthorList"]:
                last_name = author.get("LastName", "")
                fore_name = author.get("ForeName", "")
                if last_name:
                    authors.append(f"{last_name} {fore_name}".strip())

        authors_str = ", ".join(authors) if authors else "Unknown Authors"

        # Extract journal and year
        journal = article.get("Journal", {})
        journal_title = journal.get("Title", "Unknown Journal")

        year = ""
        if "JournalIssue" in journal:
            pub_date = journal["JournalIssue"].get("PubDate", {})
            year = pub_date.get("Year", "")
            if not year and "MedlineDate" in pub_date:
                import re
                year_match = re.search(r'\d{4}', pub_date["MedlineDate"])
                if year_match:
                    year = year_match.group(0)

        if not year:
            year = "Unknown Year"

        # Format the content
        content = f"""Title: {title}

Authors: {authors_str}

Journal: {journal_title}
Year: {year}
PMID: {pmid}

Abstract:
{abstract}

---
Note: Full text not available. This file contains only the abstract and metadata from PubMed.
"""

        return {
            "success": True,
            "content": content,
            "message": f"Successfully fetched abstract for PMID {pmid}"
        }

    except Exception as e:
        return {
            "success": False,
            "content": "",
            "message": f"Error fetching PubMed abstract: {str(e)}"
        }


async def download_with_retry(
    url: str,
    headers: Dict[str, str],
    max_retries: int = 5,
    initial_delay: float = 1.0
) -> requests.Response:
    """
    Download with exponential backoff retry logic.

    Args:
        url: URL to download
        headers: Request headers
        max_retries: Maximum number of retry attempts
        initial_delay: Initial delay in seconds

    Returns:
        Response object

    Raises:
        requests.exceptions.RequestException: If all retries fail
    """
    delay = initial_delay
    last_exception = None

    for attempt in range(max_retries):
        try:
            response = await asyncio.to_thread(
                requests.get,
                url,
                headers=headers,
                timeout=120
            )

            # Check for rate limiting or temporary errors
            if response.status_code in [429, 503]:
                # Check for Retry-After header
                retry_after = response.headers.get('Retry-After')
                if retry_after:
                    try:
                        delay = float(retry_after)
                    except ValueError:
                        pass

                if attempt < max_retries - 1:
                    print(f"Rate limited (attempt {attempt + 1}/{max_retries}), waiting {delay:.1f}s...")
                    await asyncio.sleep(delay)
                    delay *= 2  # Exponential backoff
                    continue

            # Success or non-retryable error
            return response

        except requests.exceptions.RequestException as e:
            last_exception = e
            if attempt < max_retries - 1:
                print(f"Request failed (attempt {attempt + 1}/{max_retries}): {str(e)}")
                await asyncio.sleep(delay)
                delay *= 2
            else:
                raise

    # Should not reach here, but just in case
    if last_exception:
        raise last_exception
    else:
        raise requests.exceptions.RequestException("Download failed after all retries")


async def download_paper_from_doi(
    doi: Optional[str] = None,
    output_path: Optional[str] = None,
    jina_api_key: Optional[str] = None,
    pmid: Optional[str] = None,
    entrez_email: Optional[str] = None,
    ncbi_api_key: Optional[str] = None,
    prefer_pmc: bool = True,
    force_download: bool = False
) -> Dict[str, Any]:
    """
    Download paper from DOI or PMC using NCBI efetch and Jina Reader API.

    At least one of DOI or PMID must be provided.

    Implements improved rate limiting, retry logic, and NCBI-compliant headers
    to prevent DDoS detection and blocking.

    Download strategy:
    1. Check cache first
    2. If PMC available: Try NCBI efetch (official API, no blocking)
    3. Fallback to DOI via Jina Reader API
    4. Final fallback to PubMed abstract

    Args:
        doi: Paper DOI (optional if pmid provided)
        output_path: Optional output file path (default: ./papers/{doi_sanitized}.txt or ./papers/pmid_{pmid}.txt)
        jina_api_key: Optional Jina API key for faster downloads
        pmid: Optional PubMed ID (optional if doi provided)
        entrez_email: Email for Entrez API (required if pmid is used for PMC check)
        ncbi_api_key: Optional NCBI API key for higher rate limits (10 req/s vs 3 req/s)
        prefer_pmc: If True, tries PMC before DOI (default: True to prioritize PMC)
        force_download: If True, re-downloads even if file exists (default: False)

    Returns:
        Dictionary with success status, file path, message, and source (PMC, DOI, or Cache)
    """
    # Validation: At least one identifier must be provided
    if not doi and not pmid:
        return {
            "success": False,
            "file_path": "",
            "source": "",
            "message": "Error: At least one of DOI or PMID must be provided"
        }

    try:
        # Check cache first unless force_download is True
        if not force_download:
            # Try checking with both DOI and PMID
            if doi:
                cached_path = check_cached_paper(doi, output_path)
                if cached_path:
                    return {
                        "success": True,
                        "file_path": cached_path,
                        "source": "Cache",
                        "message": f"Paper already downloaded at {cached_path}"
                    }
            if pmid and not doi:
                # Also check PMID-based filename (both full text and abstract)
                sanitized_pmid = pmid.replace("/", "_").replace(":", "_")
                pmid_path = f"./papers/pmid_{sanitized_pmid}.txt"
                pmid_abstract_path = f"./papers/pmid_{sanitized_pmid}_abstract.txt"

                if os.path.exists(pmid_path):
                    return {
                        "success": True,
                        "file_path": pmid_path,
                        "source": "Cache",
                        "message": f"Paper already downloaded at {pmid_path}"
                    }
                elif os.path.exists(pmid_abstract_path):
                    return {
                        "success": True,
                        "file_path": pmid_abstract_path,
                        "source": "Cache",
                        "message": f"Abstract already downloaded at {pmid_abstract_path}"
                    }

        actual_url = None
        source = "DOI"
        pmcid = None

        # Try PMC first if requested
        if prefer_pmc and pmid and entrez_email:
            # NCBI rate limiting: 2-4 seconds with random jitter
            delay = random.uniform(2.0, 4.0)
            print(f"Waiting {delay:.1f}s for NCBI rate limiting...")
            await asyncio.sleep(delay)

            pmcid = await check_pmc_availability(pmid, entrez_email, ncbi_api_key)
            if pmcid:
                print(f"Found PMC version: {pmcid}, trying NCBI efetch...")

                # Try NCBI efetch first (official API, no DDoS blocking)
                efetch_result = await fetch_pmc_fulltext(pmcid, entrez_email, ncbi_api_key)

                if efetch_result["success"]:
                    # Determine output path
                    final_output_path = output_path
                    if final_output_path is None:
                        if doi:
                            sanitized_doi = doi.replace("/", "_").replace(":", "_")
                            final_output_path = f"./papers/{sanitized_doi}.txt"
                        elif pmid:
                            sanitized_pmid = pmid.replace("/", "_").replace(":", "_")
                            final_output_path = f"./papers/pmid_{sanitized_pmid}.txt"
                        else:
                            sanitized_pmcid = pmcid.replace("/", "_").replace(":", "_")
                            final_output_path = f"./papers/{sanitized_pmcid}.txt"

                    # Create directory if needed
                    dir_path = os.path.dirname(final_output_path)
                    if dir_path:
                        os.makedirs(dir_path, exist_ok=True)

                    # Save to file
                    with open(final_output_path, "w", encoding="utf-8") as f:
                        f.write(efetch_result["content"])

                    return {
                        "success": True,
                        "file_path": final_output_path,
                        "source": "PMC_efetch",
                        "message": f"Successfully downloaded paper from PMC via NCBI efetch to {final_output_path}"
                    }
                else:
                    # efetch failed, fall back to Jina Reader
                    print(f"NCBI efetch failed: {efetch_result['message']}, falling back to Jina Reader...")
                    actual_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                    source = "PMC"

        # Use DOI by default or as fallback
        if not actual_url and doi:
            doi_url = f"https://doi.org/{doi}"
            # Add delay before DOI request too
            await asyncio.sleep(random.uniform(1.0, 2.0))
            response = await asyncio.to_thread(requests.get, doi_url, allow_redirects=True, timeout=30)
            actual_url = response.url
            print(f"Using DOI URL: {actual_url}")

        # If still no URL and only PMID provided, fetch abstract as fallback
        if not actual_url:
            if not pmid:
                return {
                    "success": False,
                    "file_path": "",
                    "source": "Unknown",
                    "message": "No DOI or PMID provided"
                }

            if not entrez_email:
                return {
                    "success": False,
                    "file_path": "",
                    "source": "PubMed_Abstract",
                    "message": "Email required for PubMed abstract fallback"
                }

            print(f"No PMC or DOI available for PMID {pmid}, fetching abstract from PubMed...")

            # Fetch abstract from PubMed
            abstract_result = await fetch_pubmed_abstract(pmid, entrez_email)

            if not abstract_result["success"]:
                return {
                    "success": False,
                    "file_path": "",
                    "source": "PubMed_Abstract",
                    "message": abstract_result["message"]
                }

            # Determine output path for abstract
            final_output_path = output_path
            if final_output_path is None:
                sanitized_pmid = pmid.replace("/", "_").replace(":", "_")
                final_output_path = f"./papers/pmid_{sanitized_pmid}_abstract.txt"

            # Create directory if needed
            dir_path = os.path.dirname(final_output_path)
            if dir_path:
                os.makedirs(dir_path, exist_ok=True)

            # Save abstract to file
            with open(final_output_path, "w", encoding="utf-8") as f:
                f.write(abstract_result["content"])

            return {
                "success": True,
                "file_path": final_output_path,
                "source": "PubMed_Abstract",
                "message": f"Successfully saved abstract to {final_output_path}. Note: Full text not available for PMID {pmid}"
            }

        # Call Jina Reader with NCBI-compliant headers
        jina_url = f"https://r.jina.ai/{actual_url}"

        # NCBI-compliant User-Agent with email (recommended by NCBI)
        user_agent = f"PubMedDownloadMCP/1.0 ({entrez_email or 'unknown@example.com'})"

        headers = {
            "User-Agent": user_agent,
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5"
        }

        # Add Jina-specific headers if API key is provided
        if jina_api_key:
            headers["Authorization"] = f"Bearer {jina_api_key}"
            headers["X-Respond-With"] = "markdown"
            headers["X-With-Generated-Alt"] = "true"
            headers["X-Remove-Selector"] = "header,footer,nav,.cookie-banner,.advertisement"

        # Download with retry logic
        try:
            result = await download_with_retry(jina_url, headers, max_retries=5)
        except requests.exceptions.RequestException as e:
            # Check if it's a SecurityCompromiseError from Jina
            error_msg = str(e)
            if "SecurityCompromiseError" in error_msg or "DDoS" in error_msg:
                return {
                    "success": False,
                    "file_path": "",
                    "source": source,
                    "message": f"PMC access blocked due to rate limiting. Please wait and try again later, or use DOI directly. Error: {error_msg}"
                }
            raise

        # Check for errors or authentication issues
        if result.status_code in [401, 403, 405]:
            # Try with minimal headers (browser-like but simpler)
            headers_basic = {
                "User-Agent": user_agent,
                "Accept": "*/*"
            }
            result = await download_with_retry(jina_url, headers_basic, max_retries=3)

        # If still getting errors, try search mode
        if result.status_code >= 400:
            try:
                content = result.text.lower()
                if any(error in content for error in ["cookies disabled", "cookie", "access denied", "403 forbidden", "paywall"]):
                    # Try search mode which might bypass restrictions
                    jina_url_search = f"https://s.jina.ai/{actual_url}"
                    headers_simple = {
                        "User-Agent": user_agent
                    }
                    result = await download_with_retry(jina_url_search, headers_simple, max_retries=3)
            except:
                pass

        result.raise_for_status()

        # Determine output path
        final_output_path = output_path
        if final_output_path is None:
            # Determine filename based on available identifiers
            if doi:
                sanitized_doi = doi.replace("/", "_").replace(":", "_")
                final_output_path = f"./papers/{sanitized_doi}.txt"
            elif pmid:
                sanitized_pmid = pmid.replace("/", "_").replace(":", "_")
                final_output_path = f"./papers/pmid_{sanitized_pmid}.txt"
            else:
                final_output_path = "./papers/unknown_paper.txt"

        # Create directory if needed
        dir_path = os.path.dirname(final_output_path)
        if dir_path:  # Only create if there's a directory path
            os.makedirs(dir_path, exist_ok=True)

        # Save to file
        with open(final_output_path, "w", encoding="utf-8") as f:
            f.write(result.text)

        return {
            "success": True,
            "file_path": final_output_path,
            "source": source,
            "message": f"Successfully downloaded paper from {source} to {final_output_path}"
        }

    except requests.exceptions.RequestException as e:
        error_msg = str(e)

        # Special handling for SecurityCompromiseError
        if "SecurityCompromiseError" in error_msg or "blocked" in error_msg.lower():
            return {
                "success": False,
                "file_path": "",
                "source": source if 'source' in locals() else "Unknown",
                "message": f"Access blocked by rate limiting. Please try again later or use prefer_pmc=False. Error: {error_msg}"
            }

        return {
            "success": False,
            "file_path": "",
            "source": source if 'source' in locals() else "Unknown",
            "message": f"Network error: {error_msg}"
        }
    except OSError as e:
        return {
            "success": False,
            "file_path": "",
            "source": source if 'source' in locals() else "Unknown",
            "message": f"File system error: {str(e)}"
        }
    except Exception as e:
        return {
            "success": False,
            "file_path": "",
            "source": source if 'source' in locals() else "Unknown",
            "message": f"Error downloading paper: {str(e)}"
        }
