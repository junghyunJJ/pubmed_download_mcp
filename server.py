"""
FastMCP server for PubMed search and paper download
"""

import os
from typing import List, Dict, Any, Optional
from fastmcp import FastMCP
from utils.pubmed_search import search_lr, search_general, search_minimal
from utils.paper_download import download_paper_from_doi


# Get configuration from environment variables
entrez_email = os.getenv("PUBMED_EMAIL", "")
jina_api_key = os.getenv("JINA_API_KEY", None)
ncbi_api_key = os.getenv("NCBI_API_KEY", None)

# Create FastMCP server
mcp = FastMCP("PubMed Download")


@mcp.tool()
async def search_pubmed_lr(
    ligand: str,
    receptor: str,
    max_results: int = 5,
    sort_by: str = "pub_date"
) -> List[Dict[str, Any]]:
    """
    Search PubMed for ligand-receptor interaction papers.

    Returns papers discussing the interaction, signaling, binding, or pathway
    between the specified ligand and receptor genes.

    Args:
        ligand: Ligand gene name (e.g., 'IL2', 'TGFB1', 'VEGFA')
        receptor: Receptor gene name (e.g., 'IL2RA', 'TGFBR1', 'KDR')
        max_results: Maximum number of results to return (default: 5)
        sort_by: Sort order - 'pub_date' for newest first, 'relevance' for most relevant (default: 'pub_date')

    Returns:
        List of papers with pmid, title, abstract, authors, year, journal, doi, pmcid, ligand, receptor
    """
    global entrez_email
    if not entrez_email:
        raise ValueError("PUBMED_EMAIL environment variable or --email argument required")

    return await search_lr(
        ligand=ligand,
        receptor=receptor,
        max_results=max_results,
        sort_by=sort_by,
        email=entrez_email
    )


@mcp.tool()
async def search_pubmed(
    query: str,
    max_results: int = 10,
    sort_by: str = "pub_date"
) -> List[Dict[str, Any]]:
    """
    Search PubMed with a general query.

    Supports standard PubMed query syntax including field tags, Boolean operators, and filters.

    Args:
        query: PubMed search query (e.g., 'cancer AND immunotherapy', 'BRCA1[Gene] AND breast cancer')
        max_results: Maximum number of results to return (default: 10)
        sort_by: Sort order - 'pub_date' for newest first, 'relevance' for most relevant (default: 'pub_date')

    Returns:
        List of papers with pmid, title, abstract, authors, year, journal, doi, pmcid
    """
    global entrez_email
    if not entrez_email:
        raise ValueError("PUBMED_EMAIL environment variable or --email argument required")

    return await search_general(
        query=query,
        max_results=max_results,
        sort_by=sort_by,
        email=entrez_email
    )


@mcp.tool()
async def search_pubmed_minimal(
    query: str,
    max_results: int = 10,
    sort_by: str = "pub_date"
) -> List[Dict[str, Any]]:
    """
    Search PubMed with minimal metadata (no abstract) for token efficiency.

    Returns only essential information: pmid, title, first_author, year, journal, doi, pmcid.
    This reduces token usage by ~70-80% compared to search_pubmed.

    Supports standard PubMed query syntax including field tags, Boolean operators, and filters.

    Args:
        query: PubMed search query (e.g., 'cancer AND immunotherapy', 'BRCA1[Gene] AND breast cancer')
        max_results: Maximum number of results to return (default: 10)
        sort_by: Sort order - 'pub_date' for newest first, 'relevance' for most relevant (default: 'pub_date')

    Returns:
        List of papers with pmid, title, first_author, year, journal, doi, pmcid (no abstract)
    """
    global entrez_email
    if not entrez_email:
        raise ValueError("PUBMED_EMAIL environment variable or --email argument required")

    return await search_minimal(
        query=query,
        max_results=max_results,
        sort_by=sort_by,
        email=entrez_email
    )


@mcp.tool()
async def download_paper(
    doi: Optional[str] = None,
    output_path: Optional[str] = None,
    pmid: Optional[str] = None,
    prefer_pmc: bool = True
) -> Dict[str, Any]:
    """
    Download paper content from PubMed Central (PMC) or DOI using Jina Reader API.

    At least one of DOI or PMID must be provided.

    By default prioritizes PMC (prefer_pmc=True) for better content extraction.
    Implements caching, exponential backoff retry, and NCBI-compliant headers to prevent DDoS detection.

    Args:
        doi: Paper DOI (e.g., '10.1038/nature12345') - optional if pmid provided
        output_path: Optional output file path (default: ./papers/{doi_sanitized}.txt or ./papers/pmid_{pmid}.txt)
        pmid: Optional PubMed ID - optional if doi provided
        prefer_pmc: If True, tries PMC before DOI (default: True to prioritize PMC downloads)

    Returns:
        Dictionary with success status, file path, source (PMC/DOI/Cache), and message
    """
    global entrez_email, jina_api_key, ncbi_api_key

    return await download_paper_from_doi(
        doi=doi,
        output_path=output_path,
        jina_api_key=jina_api_key,
        pmid=pmid,
        entrez_email=entrez_email,
        ncbi_api_key=ncbi_api_key,
        prefer_pmc=prefer_pmc
    )


def set_config(email: str, jina_key: Optional[str] = None, ncbi_key: Optional[str] = None):
    """Set configuration for the server"""
    global entrez_email, jina_api_key, ncbi_api_key
    entrez_email = email
    if jina_key:
        jina_api_key = jina_key
    if ncbi_key:
        ncbi_api_key = ncbi_key


if __name__ == "__main__":
    mcp.run()
