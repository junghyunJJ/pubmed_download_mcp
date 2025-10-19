"""
FastMCP server for PubMed search and paper download
"""

import os
import argparse
from typing import List, Dict, Any, Optional
from fastmcp import FastMCP
from utils.pubmed_search import search_lr, search_general
from utils.paper_download import download_paper_from_doi


# Parse command-line arguments
parser = argparse.ArgumentParser(description="PubMed Download MCP Server")
parser.add_argument("--email", type=str, help="Email address for NCBI Entrez API (required)")
parser.add_argument("--jina-api-key", type=str, help="Jina AI API key for paper downloads (optional)")

# Get configuration from command-line arguments or environment variables
args, _ = parser.parse_known_args()
entrez_email = args.email or os.getenv("PUBMED_EMAIL", "")
jina_api_key = args.jina_api_key or os.getenv("JINA_API_KEY", None)

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
async def download_paper(
    doi: str,
    output_path: Optional[str] = None,
    pmid: Optional[str] = None
) -> Dict[str, Any]:
    """
    Download paper content from PubMed Central (PMC) or DOI using Jina Reader API.

    Tries PMC first if PMID is provided (more reliable), falls back to DOI if PMC unavailable.
    Saves the paper as a text file with full content including title, authors, abstract, and body text.

    Args:
        doi: Paper DOI (e.g., '10.1038/nature12345')
        output_path: Optional output file path (default: ./papers/{doi_sanitized}.txt)
        pmid: Optional PubMed ID to check PMC availability (recommended for better success rate)

    Returns:
        Dictionary with success status, file path, source (PMC or DOI), and message
    """
    global entrez_email, jina_api_key

    return await download_paper_from_doi(
        doi=doi,
        output_path=output_path,
        jina_api_key=jina_api_key,
        pmid=pmid,
        entrez_email=entrez_email
    )


def set_config(email: str, jina_key: Optional[str] = None):
    """Set configuration for the server"""
    global entrez_email, jina_api_key
    entrez_email = email
    if jina_key:
        jina_api_key = jina_key


if __name__ == "__main__":
    mcp.run()
