"""
⚠️ DEPRECATED: This module is no longer used as of v2.2

PMC changed their website to use JavaScript rendering, making HTML parsing unreliable.
The Jina Reader API now handles PMC content extraction directly in paper_download.py.

This file is kept for reference but is not imported or used by the server.

---

PMC (PubMed Central) HTML parser for extracting paper content.
Converts PMC HTML articles to clean text format with NCBI-compliant headers and rate limiting.
"""

import asyncio
import random
import time
from typing import Optional, Dict, Any
import requests
from bs4 import BeautifulSoup


async def download_pmc_article_async(
    pmc_id: str,
    email: str,
    api_key: Optional[str] = None,
    include_references: bool = False,
    max_retries: int = 5
) -> Dict[str, Any]:
    """
    Download and parse a PMC article to clean text format.

    Args:
        pmc_id: PMC ID (with or without 'PMC' prefix)
        email: Email for NCBI compliance (required)
        api_key: Optional NCBI API key for higher rate limits
        include_references: Whether to include references section in output
        max_retries: Maximum number of retry attempts

    Returns:
        Dict with keys:
            - success: bool indicating success
            - content: parsed text content if successful
            - message: status or error message
            - source: "pmc_html"
    """
    # Clean PMC ID (remove 'PMC' prefix if present)
    if pmc_id.startswith('PMC'):
        pmc_id = pmc_id[3:]

    # Wrap synchronous function in async
    return await asyncio.to_thread(
        download_pmc_article_with_retry,
        pmc_id,
        email,
        api_key,
        include_references,
        max_retries
    )


def download_pmc_article_with_retry(
    pmc_id: str,
    email: str,
    api_key: Optional[str] = None,
    include_references: bool = False,
    max_retries: int = 5
) -> Dict[str, Any]:
    """
    Download PMC article with retry logic and NCBI-compliant headers.
    """
    url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/"

    # NCBI-compliant headers
    headers = {
        'User-Agent': f'PubMedDownloadMCP/1.0 ({email})',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Connection': 'keep-alive',
    }

    # Add API key header if provided (for higher rate limits)
    if api_key:
        headers['api_key'] = api_key

    for attempt in range(max_retries):
        try:
            # Rate limiting: 2-4 second random delay to prevent DDoS detection
            if attempt > 0:
                # Exponential backoff on retries
                delay = min(2 ** attempt, 16)  # Cap at 16 seconds
                time.sleep(delay)
            else:
                # Initial request delay
                delay = random.uniform(2, 4)
                time.sleep(delay)

            # Make request with NCBI-compliant headers
            response = requests.get(url, headers=headers, timeout=30)

            # Handle rate limiting
            if response.status_code == 429:
                # Check for Retry-After header
                retry_after = response.headers.get('Retry-After', str(5 * (attempt + 1)))
                try:
                    wait_time = int(retry_after)
                except ValueError:
                    wait_time = 5 * (attempt + 1)

                if attempt < max_retries - 1:
                    time.sleep(wait_time)
                    continue
                else:
                    return {
                        'success': False,
                        'message': f'Rate limited after {max_retries} attempts',
                        'source': 'pmc_html'
                    }

            # Handle server errors with retry
            if response.status_code in [500, 502, 503, 504]:
                if attempt < max_retries - 1:
                    continue
                else:
                    return {
                        'success': False,
                        'message': f'Server error {response.status_code} after {max_retries} attempts',
                        'source': 'pmc_html'
                    }

            # Check for successful response
            response.raise_for_status()

            # Parse the HTML content
            parsed_content = parse_pmc_html(response.text, include_references)

            if parsed_content:
                return {
                    'success': True,
                    'content': parsed_content,
                    'message': f'Successfully downloaded PMC{pmc_id} via HTML parser',
                    'source': 'pmc_html'
                }
            else:
                return {
                    'success': False,
                    'message': 'Failed to parse PMC HTML content',
                    'source': 'pmc_html'
                }

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 403:
                return {
                    'success': False,
                    'message': f'Access forbidden (403) for PMC{pmc_id}. Paper may be embargoed or restricted.',
                    'source': 'pmc_html'
                }
            elif e.response.status_code == 404:
                return {
                    'success': False,
                    'message': f'PMC{pmc_id} not found (404)',
                    'source': 'pmc_html'
                }
            elif attempt < max_retries - 1:
                continue
            else:
                return {
                    'success': False,
                    'message': f'HTTP error {e.response.status_code}: {str(e)}',
                    'source': 'pmc_html'
                }

        except requests.exceptions.Timeout:
            if attempt < max_retries - 1:
                continue
            else:
                return {
                    'success': False,
                    'message': 'Request timeout after multiple attempts',
                    'source': 'pmc_html'
                }

        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                continue
            else:
                return {
                    'success': False,
                    'message': f'Request error: {str(e)}',
                    'source': 'pmc_html'
                }
        except Exception as e:
            return {
                'success': False,
                'message': f'Unexpected error: {str(e)}',
                'source': 'pmc_html'
            }

    return {
        'success': False,
        'message': f'Failed to download PMC{pmc_id} after {max_retries} attempts',
        'source': 'pmc_html'
    }


def parse_pmc_html(html_content: str, include_references: bool = False) -> Optional[str]:
    """
    Parse PMC HTML content and extract clean text.

    Args:
        html_content: Raw HTML from PMC
        include_references: Whether to include references section

    Returns:
        Formatted text content or None if parsing fails
    """
    try:
        soup = BeautifulSoup(html_content, 'html.parser')

        # Extract title
        title = ""
        title_elem = soup.find('h1', class_='content-title')
        if not title_elem:
            # Try alternative selectors
            title_elem = soup.find('h1', {'id': 'article-title'})
            if not title_elem:
                title_elem = soup.find('div', class_='fm-title')
        if title_elem:
            title = title_elem.get_text(strip=True)

        # Extract authors
        authors = []
        author_list = soup.find('div', class_='contrib-group')
        if not author_list:
            author_list = soup.find('div', class_='fm-author')
        if author_list:
            for author in author_list.find_all(['span', 'a'], class_=['contrib-name', 'contrib']):
                author_text = author.get_text(strip=True)
                if author_text:
                    authors.append(author_text)

        # Extract abstract
        abstract = ""
        abstract_elem = soup.find('div', class_='abstract')
        if not abstract_elem:
            abstract_elem = soup.find('div', {'id': 'abstract'})
        if abstract_elem:
            # Remove any sub-headers within abstract
            for header in abstract_elem.find_all(['h2', 'h3', 'h4']):
                header.decompose()
            abstract = abstract_elem.get_text(separator=' ', strip=True)

        # Extract main body
        body_sections = []

        # Look for main content area
        main_content = soup.find('div', class_='body')
        if not main_content:
            main_content = soup.find('div', {'id': 'body'})
        if not main_content:
            main_content = soup.find('main')
        if not main_content:
            # Try to find sections directly
            main_content = soup

        # Find all sections in the body
        for section in main_content.find_all(['div', 'section'], class_=['sec', 'section']):
            # Skip references section if not requested
            section_id = section.get('id', '').lower()
            section_class = ' '.join(section.get('class', [])).lower()

            if not include_references and ('ref' in section_id or 'reference' in section_id or
                                          'ref' in section_class or 'reference' in section_class):
                continue

            # Get section title if available
            section_title = None
            title_elem = section.find(['h2', 'h3', 'h4'])
            if title_elem:
                section_title = title_elem.get_text(strip=True)
                title_elem.decompose()  # Remove to avoid duplication

            # Get section content
            section_text = section.get_text(separator=' ', strip=True)

            if section_text:
                if section_title:
                    body_sections.append(f"\n## {section_title}\n\n{section_text}")
                else:
                    body_sections.append(section_text)

        # Extract references if requested
        references = ""
        if include_references:
            ref_section = soup.find('div', {'id': 'references'})
            if not ref_section:
                ref_section = soup.find('div', class_='ref-list')
            if not ref_section:
                ref_section = soup.find('section', {'id': 'references'})

            if ref_section:
                ref_items = []
                for ref in ref_section.find_all(['li', 'div'], class_=['ref', 'ref-cit']):
                    ref_text = ref.get_text(separator=' ', strip=True)
                    if ref_text:
                        ref_items.append(ref_text)

                if ref_items:
                    references = "\n\nReferences:\n" + "\n\n".join(ref_items)

        # Format output
        output_parts = []

        if title:
            output_parts.append(f"Title: {title}")

        if authors:
            output_parts.append(f"\nAuthors: {', '.join(authors)}")

        if abstract:
            output_parts.append(f"\n\nAbstract:\n{abstract}")

        if body_sections:
            output_parts.append(f"\n\nBody:\n{''.join(body_sections)}")

        if references:
            output_parts.append(references)

        # Return formatted content if we found meaningful content
        if output_parts:
            return '\n'.join(output_parts)
        else:
            return None

    except Exception as e:
        print(f"Error parsing PMC HTML: {str(e)}")
        return None


# Test function
if __name__ == "__main__":
    import os

    # Test with PMC12165575 (PMID: 39945774)
    email = os.getenv("PUBMED_EMAIL", "test@example.com")
    api_key = os.getenv("NCBI_API_KEY")

    async def test():
        print("Testing PMC download with NCBI-compliant headers...")
        result = await download_pmc_article_async(
            "12165575",
            email,
            api_key,
            include_references=True
        )

        if result['success']:
            print("Success! Content preview:")
            print(result['content'][:500])

            # Save to file
            with open("PMC12165575_test.txt", "w") as f:
                f.write(result['content'])
            print("\nFull content saved to PMC12165575_test.txt")
        else:
            print(f"Failed: {result['message']}")

    asyncio.run(test())
