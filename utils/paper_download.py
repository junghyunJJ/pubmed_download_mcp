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


async def _try_jina_download(
    url: str,
    jina_api_key: Optional[str],
    entrez_email: Optional[str]
) -> Optional[str]:
    """
    Helper function to download content via Jina Reader API.

    Args:
        url: Target URL to fetch
        jina_api_key: Optional Jina API key
        entrez_email: Email for User-Agent

    Returns:
        Content string if successful, None if failed
    """
    jina_url = f"https://r.jina.ai/{url}"
    user_agent = f"PubMedDownloadMCP/1.0 ({entrez_email or 'unknown@example.com'})"

    headers = {
        "User-Agent": user_agent,
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.5"
    }

    if jina_api_key:
        headers["Authorization"] = f"Bearer {jina_api_key}"
        headers["X-Respond-With"] = "markdown"
        headers["X-With-Generated-Alt"] = "true"
        headers["X-Remove-Selector"] = "header,footer,nav,.cookie-banner,.advertisement"

    try:
        result = await download_with_retry(jina_url, headers, max_retries=5)

        # Check for errors or authentication issues
        if result.status_code in [401, 403, 405]:
            headers_basic = {"User-Agent": user_agent, "Accept": "*/*"}
            result = await download_with_retry(jina_url, headers_basic, max_retries=3)

        # If still getting errors, try search mode
        if result.status_code >= 400:
            try:
                content = result.text.lower()
                if any(error in content for error in ["cookies disabled", "cookie", "access denied", "403 forbidden", "paywall"]):
                    jina_url_search = f"https://s.jina.ai/{url}"
                    headers_simple = {"User-Agent": user_agent}
                    result = await download_with_retry(jina_url_search, headers_simple, max_retries=3)
            except:
                pass

        result.raise_for_status()
        return result.text

    except requests.exceptions.RequestException as e:
        print(f"Jina download failed for {url}: {str(e)}")
        return None


def _determine_output_path(
    output_path: Optional[str],
    doi: Optional[str],
    pmid: Optional[str],
    pmcid: Optional[str] = None,
    is_abstract: bool = False
) -> str:
    """Helper function to determine the output file path."""
    if output_path:
        return output_path

    if doi:
        sanitized_doi = doi.replace("/", "_").replace(":", "_")
        return f"./papers/{sanitized_doi}.txt"
    elif pmid:
        sanitized_pmid = pmid.replace("/", "_").replace(":", "_")
        suffix = "_abstract" if is_abstract else ""
        return f"./papers/pmid_{sanitized_pmid}{suffix}.txt"
    elif pmcid:
        sanitized_pmcid = pmcid.replace("/", "_").replace(":", "_")
        return f"./papers/{sanitized_pmcid}.txt"
    else:
        return "./papers/unknown_paper.txt"


def _save_content(content: str, file_path: str) -> None:
    """Helper function to save content to file."""
    dir_path = os.path.dirname(file_path)
    if dir_path:
        os.makedirs(dir_path, exist_ok=True)
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(content)


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

    Download strategy (v4.1 - improved fallback):
    1. Check cache first
    2. If PMC available: Try NCBI efetch → If fails, try PMC via Jina
    3. If PMC fails or unavailable: Try DOI via Jina Reader API
    4. Final fallback: PubMed abstract only

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
        Dictionary with success status, file path, message, and source (PMC_efetch/PMC/DOI/PubMed_Abstract/Cache)
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
        # ===== STEP 1: Check cache =====
        if not force_download:
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

        pmcid = None

        # ===== STEP 2: Try PMC if preferred =====
        if prefer_pmc and pmid and entrez_email:
            # NCBI rate limiting
            delay = random.uniform(2.0, 4.0)
            print(f"Waiting {delay:.1f}s for NCBI rate limiting...")
            await asyncio.sleep(delay)

            pmcid = await check_pmc_availability(pmid, entrez_email, ncbi_api_key)
            if pmcid:
                print(f"Found PMC version: {pmcid}, trying NCBI efetch...")

                # 2a. Try NCBI efetch first (official API, no DDoS blocking)
                efetch_result = await fetch_pmc_fulltext(pmcid, entrez_email, ncbi_api_key)

                if efetch_result["success"]:
                    final_path = _determine_output_path(output_path, doi, pmid, pmcid)
                    _save_content(efetch_result["content"], final_path)
                    return {
                        "success": True,
                        "file_path": final_path,
                        "source": "PMC_efetch",
                        "message": f"Successfully downloaded paper from PMC via NCBI efetch to {final_path}"
                    }

                # 2b. efetch failed, try PMC via Jina Reader
                print(f"NCBI efetch failed: {efetch_result['message']}, trying PMC via Jina...")
                pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                pmc_content = await _try_jina_download(pmc_url, jina_api_key, entrez_email)

                if pmc_content:
                    final_path = _determine_output_path(output_path, doi, pmid, pmcid)
                    _save_content(pmc_content, final_path)
                    return {
                        "success": True,
                        "file_path": final_path,
                        "source": "PMC",
                        "message": f"Successfully downloaded paper from PMC via Jina to {final_path}"
                    }

                print(f"PMC Jina also failed, falling back to DOI...")

        # ===== STEP 3: Try DOI =====
        if doi:
            print(f"Trying DOI: {doi}")
            await asyncio.sleep(random.uniform(1.0, 2.0))

            try:
                response = await asyncio.to_thread(
                    requests.get, f"https://doi.org/{doi}",
                    allow_redirects=True, timeout=30
                )
                actual_url = response.url
                print(f"Resolved DOI URL: {actual_url}")

                doi_content = await _try_jina_download(actual_url, jina_api_key, entrez_email)

                if doi_content:
                    final_path = _determine_output_path(output_path, doi, pmid)
                    _save_content(doi_content, final_path)
                    return {
                        "success": True,
                        "file_path": final_path,
                        "source": "DOI",
                        "message": f"Successfully downloaded paper from DOI to {final_path}"
                    }

                print(f"DOI download failed, falling back to abstract...")

            except requests.exceptions.RequestException as e:
                print(f"DOI resolution failed: {str(e)}, falling back to abstract...")

        # ===== STEP 4: Fallback to PubMed abstract =====
        if pmid and entrez_email:
            print(f"Fetching abstract for PMID {pmid}...")
            abstract_result = await fetch_pubmed_abstract(pmid, entrez_email)

            if abstract_result["success"]:
                final_path = _determine_output_path(output_path, doi, pmid, is_abstract=True)
                _save_content(abstract_result["content"], final_path)
                return {
                    "success": True,
                    "file_path": final_path,
                    "source": "PubMed_Abstract",
                    "message": f"Full text unavailable. Saved abstract to {final_path}"
                }

        # ===== STEP 5: All methods failed =====
        methods_tried = []
        if prefer_pmc and pmid:
            methods_tried.append("PMC efetch")
            if pmcid:
                methods_tried.append("PMC Jina")
        if doi:
            methods_tried.append("DOI")
        if pmid:
            methods_tried.append("Abstract")

        return {
            "success": False,
            "file_path": "",
            "source": "",
            "message": f"All download methods failed. Tried: {', '.join(methods_tried) or 'none'}"
        }

    except requests.exceptions.RequestException as e:
        error_msg = str(e)
        if "SecurityCompromiseError" in error_msg or "blocked" in error_msg.lower():
            return {
                "success": False,
                "file_path": "",
                "source": "",
                "message": f"Access blocked by rate limiting. Please try again later. Error: {error_msg}"
            }
        return {
            "success": False,
            "file_path": "",
            "source": "",
            "message": f"Network error: {error_msg}"
        }
    except OSError as e:
        return {
            "success": False,
            "file_path": "",
            "source": "",
            "message": f"File system error: {str(e)}"
        }
    except Exception as e:
        return {
            "success": False,
            "file_path": "",
            "source": "",
            "message": f"Error downloading paper: {str(e)}"
        }
