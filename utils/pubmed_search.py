"""
PubMed search functionality using Biopython Entrez API
"""

import re
import asyncio
from typing import List, Dict, Any
from Bio import Entrez


async def search_lr(
    ligand: str,
    receptor: str,
    max_results: int,
    sort_by: str,
    email: str
) -> List[Dict[str, Any]]:
    """
    Search PubMed for ligand-receptor interaction papers.

    Args:
        ligand: Ligand gene name
        receptor: Receptor gene name
        max_results: Maximum number of results
        sort_by: Sort order ('pub_date' or 'relevance')
        email: Email for NCBI Entrez API

    Returns:
        List of paper information dictionaries
    """
    Entrez.email = email

    # Construct search query
    query = f"({ligand}[Gene Symbol] OR {ligand}[MeSH Terms]) AND ({receptor}[Gene Symbol] OR {receptor}[MeSH Terms]) AND (interaction OR signaling OR binding OR pathway)"

    try:
        # Search PubMed (wrapped in asyncio.to_thread for non-blocking)
        def _search():
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort=sort_by
            )
            search_results = Entrez.read(handle)
            handle.close()
            return search_results

        search_results = await asyncio.to_thread(_search)
        id_list = search_results["IdList"]

        if not id_list:
            return []

        # Fetch paper details (wrapped in asyncio.to_thread for non-blocking)
        await asyncio.sleep(0.5)  # Rate limiting

        def _fetch():
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(id_list),
                retmode="xml"
            )
            records = Entrez.read(handle)
            handle.close()
            return records

        records = await asyncio.to_thread(_fetch)
        pubmed_articles = records.get("PubmedArticle", [])

        results = []
        for record in pubmed_articles:
            try:
                medline_citation = record.get("MedlineCitation", {})
                pubmed_data = record.get("PubmedData", {})
                article = medline_citation.get("Article", {})

                # Extract PMID
                pmid = str(medline_citation.get("PMID", ""))

                # Extract title
                title = article.get("ArticleTitle", "")

                # Extract abstract
                abstract = ""
                if "Abstract" in article:
                    abstract_texts = article["Abstract"].get("AbstractText", [])
                    if isinstance(abstract_texts, list):
                        abstract = " ".join(str(text) for text in abstract_texts)
                    else:
                        abstract = str(abstract_texts)

                # Extract authors
                authors = []
                if "AuthorList" in article:
                    for author in article["AuthorList"]:
                        last_name = author.get("LastName", "")
                        fore_name = author.get("ForeName", "")
                        if last_name:
                            authors.append(f"{last_name} {fore_name}".strip())

                # Extract year
                year = ""
                journal = article.get("Journal", {})
                if "JournalIssue" in journal:
                    pub_date = journal["JournalIssue"].get("PubDate", {})
                    year = pub_date.get("Year", "")
                    if not year and "MedlineDate" in pub_date:
                        year_match = re.search(r'\d{4}', pub_date["MedlineDate"])
                        if year_match:
                            year = year_match.group(0)

                # Extract journal name
                journal_title = journal.get("Title", "")

                # Extract DOI and PMCID
                doi = ""
                pmcid = ""
                article_id_list = pubmed_data.get("ArticleIdList", [])
                for article_id in article_id_list:
                    id_type = article_id.attributes.get("IdType")
                    if id_type == "doi":
                        doi = str(article_id)
                    elif id_type == "pmc":
                        pmcid = str(article_id)

                article_info = {
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "authors": ", ".join(authors),
                    "year": year,
                    "journal": journal_title,
                    "doi": doi,
                    "pmcid": pmcid,
                    "ligand": ligand,
                    "receptor": receptor
                }

                results.append(article_info)

            except Exception as e:
                print(f"Error extracting paper info: {str(e)}")
                continue

        return results

    except Exception as e:
        raise Exception(f"PubMed search error: {str(e)}")


async def search_general(
    query: str,
    max_results: int,
    sort_by: str,
    email: str
) -> List[Dict[str, Any]]:
    """
    Search PubMed with a general query.

    Args:
        query: PubMed search query
        max_results: Maximum number of results
        sort_by: Sort order ('pub_date' or 'relevance')
        email: Email for NCBI Entrez API

    Returns:
        List of paper information dictionaries
    """
    Entrez.email = email

    try:
        # Search PubMed (wrapped in asyncio.to_thread for non-blocking)
        def _search():
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort=sort_by
            )
            search_results = Entrez.read(handle)
            handle.close()
            return search_results

        search_results = await asyncio.to_thread(_search)
        id_list = search_results["IdList"]

        if not id_list:
            return []

        # Fetch paper details (wrapped in asyncio.to_thread for non-blocking)
        await asyncio.sleep(0.5)  # Rate limiting

        def _fetch():
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(id_list),
                retmode="xml"
            )
            records = Entrez.read(handle)
            handle.close()
            return records

        records = await asyncio.to_thread(_fetch)
        pubmed_articles = records.get("PubmedArticle", [])

        results = []
        for record in pubmed_articles:
            try:
                medline_citation = record.get("MedlineCitation", {})
                pubmed_data = record.get("PubmedData", {})
                article = medline_citation.get("Article", {})

                # Extract PMID
                pmid = str(medline_citation.get("PMID", ""))

                # Extract title
                title = article.get("ArticleTitle", "")

                # Extract abstract
                abstract = ""
                if "Abstract" in article:
                    abstract_texts = article["Abstract"].get("AbstractText", [])
                    if isinstance(abstract_texts, list):
                        abstract = " ".join(str(text) for text in abstract_texts)
                    else:
                        abstract = str(abstract_texts)

                # Extract authors
                authors = []
                if "AuthorList" in article:
                    for author in article["AuthorList"]:
                        last_name = author.get("LastName", "")
                        fore_name = author.get("ForeName", "")
                        if last_name:
                            authors.append(f"{last_name} {fore_name}".strip())

                # Extract year
                year = ""
                journal = article.get("Journal", {})
                if "JournalIssue" in journal:
                    pub_date = journal["JournalIssue"].get("PubDate", {})
                    year = pub_date.get("Year", "")
                    if not year and "MedlineDate" in pub_date:
                        year_match = re.search(r'\d{4}', pub_date["MedlineDate"])
                        if year_match:
                            year = year_match.group(0)

                # Extract journal name
                journal_title = journal.get("Title", "")

                # Extract DOI and PMCID
                doi = ""
                pmcid = ""
                article_id_list = pubmed_data.get("ArticleIdList", [])
                for article_id in article_id_list:
                    id_type = article_id.attributes.get("IdType")
                    if id_type == "doi":
                        doi = str(article_id)
                    elif id_type == "pmc":
                        pmcid = str(article_id)

                article_info = {
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "authors": ", ".join(authors),
                    "year": year,
                    "journal": journal_title,
                    "doi": doi,
                    "pmcid": pmcid
                }

                results.append(article_info)

            except Exception as e:
                print(f"Error extracting paper info: {str(e)}")
                continue

        return results

    except Exception as e:
        raise Exception(f"PubMed search error: {str(e)}")


async def search_minimal(
    query: str,
    max_results: int,
    sort_by: str,
    email: str
) -> List[Dict[str, Any]]:
    """
    Search PubMed with minimal metadata (no abstract) for token efficiency.

    Args:
        query: PubMed search query
        max_results: Maximum number of results
        sort_by: Sort order ('pub_date' or 'relevance')
        email: Email for NCBI Entrez API

    Returns:
        List of paper information dictionaries (without abstracts)
    """
    Entrez.email = email

    try:
        # Search PubMed (wrapped in asyncio.to_thread for non-blocking)
        def _search():
            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=max_results,
                sort=sort_by
            )
            search_results = Entrez.read(handle)
            handle.close()
            return search_results

        search_results = await asyncio.to_thread(_search)
        id_list = search_results["IdList"]

        if not id_list:
            return []

        # Fetch paper details (wrapped in asyncio.to_thread for non-blocking)
        await asyncio.sleep(0.5)  # Rate limiting

        def _fetch():
            handle = Entrez.efetch(
                db="pubmed",
                id=",".join(id_list),
                retmode="xml"
            )
            records = Entrez.read(handle)
            handle.close()
            return records

        records = await asyncio.to_thread(_fetch)
        pubmed_articles = records.get("PubmedArticle", [])

        results = []
        for record in pubmed_articles:
            try:
                medline_citation = record.get("MedlineCitation", {})
                pubmed_data = record.get("PubmedData", {})
                article = medline_citation.get("Article", {})

                # Extract PMID
                pmid = str(medline_citation.get("PMID", ""))

                # Extract title
                title = article.get("ArticleTitle", "")

                # Extract first author only for efficiency
                first_author = ""
                if "AuthorList" in article and len(article["AuthorList"]) > 0:
                    author = article["AuthorList"][0]
                    last_name = author.get("LastName", "")
                    fore_name = author.get("ForeName", "")
                    if last_name:
                        first_author = f"{last_name} {fore_name}".strip()

                # Extract year
                year = ""
                journal = article.get("Journal", {})
                if "JournalIssue" in journal:
                    pub_date = journal["JournalIssue"].get("PubDate", {})
                    year = pub_date.get("Year", "")
                    if not year and "MedlineDate" in pub_date:
                        year_match = re.search(r'\d{4}', pub_date["MedlineDate"])
                        if year_match:
                            year = year_match.group(0)

                # Extract journal name
                journal_title = journal.get("Title", "")

                # Extract DOI and PMCID
                doi = ""
                pmcid = ""
                article_id_list = pubmed_data.get("ArticleIdList", [])
                for article_id in article_id_list:
                    id_type = article_id.attributes.get("IdType")
                    if id_type == "doi":
                        doi = str(article_id)
                    elif id_type == "pmc":
                        pmcid = str(article_id)

                article_info = {
                    "pmid": pmid,
                    "title": title,
                    "first_author": first_author,
                    "year": year,
                    "journal": journal_title,
                    "doi": doi,
                    "pmcid": pmcid
                }

                results.append(article_info)

            except Exception as e:
                print(f"Error extracting paper info: {str(e)}")
                continue

        return results

    except Exception as e:
        raise Exception(f"PubMed search error: {str(e)}")
