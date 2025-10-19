"""
Paper download functionality using Jina Reader API
"""

import os
import asyncio
from typing import Dict, Any, Optional
import requests
from Bio import Entrez


async def check_pmc_availability(pmid: str, email: str) -> Optional[str]:
    """
    Check if a paper is available in PubMed Central.

    Args:
        pmid: PubMed ID
        email: Email for Entrez API

    Returns:
        PMC ID (e.g., "PMC1234567") if available, None otherwise
    """
    try:
        Entrez.email = email

        def _check():
            handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
            result = Entrez.read(handle)
            handle.close()
            return result

        result = await asyncio.to_thread(_check)

        # Check if PMC link exists
        if result and result[0].get("LinkSetDb"):
            for linkset in result[0]["LinkSetDb"]:
                if linkset["DbTo"] == "pmc":
                    pmc_ids = linkset.get("Link", [])
                    if pmc_ids:
                        return f"PMC{pmc_ids[0]['Id']}"
        return None
    except Exception as e:
        print(f"Error checking PMC availability: {str(e)}")
        return None


async def download_paper_from_doi(
    doi: str,
    output_path: Optional[str] = None,
    jina_api_key: Optional[str] = None,
    pmid: Optional[str] = None,
    entrez_email: Optional[str] = None
) -> Dict[str, Any]:
    """
    Download paper from DOI or PMC using Jina Reader API.

    Tries PMC first if PMID is provided, falls back to DOI if PMC is not available.

    Args:
        doi: Paper DOI
        output_path: Optional output file path (default: ./papers/{doi_sanitized}.txt)
        jina_api_key: Optional Jina API key for faster downloads
        pmid: Optional PubMed ID for checking PMC availability
        entrez_email: Email for Entrez API (required if pmid is provided)

    Returns:
        Dictionary with success status, file path, message, and source (PMC or DOI)
    """
    try:
        actual_url = None
        source = "DOI"

        # Try PMC first if PMID is provided
        if pmid and entrez_email:
            await asyncio.sleep(0.5)  # Rate limiting
            pmcid = await check_pmc_availability(pmid, entrez_email)
            if pmcid:
                actual_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
                source = "PMC"
                print(f"Found PMC version: {pmcid}")

        # Fall back to DOI if PMC not available
        if not actual_url:
            doi_url = f"https://doi.org/{doi}"
            response = await asyncio.to_thread(requests.get, doi_url, allow_redirects=True, timeout=30)
            actual_url = response.url
            print(f"Using DOI URL: {actual_url}")

        # Call Jina Reader with options to bypass cookie walls
        # Add X-Respond-With header for cleaner markdown output
        # Add X-With-Generated-Alt for better image descriptions
        jina_url = f"https://r.jina.ai/{actual_url}"

        headers = {
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5"
        }

        # Add Jina-specific headers if API key is provided
        if jina_api_key:
            headers["Authorization"] = f"Bearer {jina_api_key}"
            headers["X-Respond-With"] = "markdown"
            headers["X-With-Generated-Alt"] = "true"
            headers["X-Remove-Selector"] = "header,footer,nav,.cookie-banner,.advertisement"

        result = await asyncio.to_thread(requests.get, jina_url, headers=headers, timeout=120)

        # Check for errors or authentication issues
        if result.status_code in [401, 403, 405]:
            # Try with minimal headers (browser-like but simpler)
            headers_basic = {
                "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
                "Accept": "*/*"
            }
            result = await asyncio.to_thread(requests.get, jina_url, headers=headers_basic, timeout=120)

        # If still getting errors, try without raising exception first
        if result.status_code >= 400:
            # Check if we can read the content to see if it's a paywall/cookie issue
            try:
                content = result.text.lower()
                if any(error in content for error in ["cookies disabled", "cookie", "access denied", "403 forbidden", "paywall"]):
                    # Try search mode which might bypass restrictions
                    jina_url_search = f"https://s.jina.ai/{actual_url}"
                    headers_simple = {
                        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36"
                    }
                    result = await asyncio.to_thread(requests.get, jina_url_search, headers=headers_simple, timeout=120)
            except:
                pass

        result.raise_for_status()

        # Determine output path
        final_output_path = output_path
        if final_output_path is None:
            # Sanitize DOI for filename
            sanitized_doi = doi.replace("/", "_").replace(":", "_")
            final_output_path = f"./papers/{sanitized_doi}.txt"

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
        return {
            "success": False,
            "file_path": "",
            "message": f"Network error: {str(e)}"
        }
    except OSError as e:
        return {
            "success": False,
            "file_path": "",
            "message": f"File system error: {str(e)}"
        }
    except Exception as e:
        return {
            "success": False,
            "file_path": "",
            "message": f"Error downloading paper: {str(e)}"
        }
