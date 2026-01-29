# Methods Section: PubMed Download MCP Server

## Development of PubMed Literature Retrieval System

### Overview

We developed a Model Context Protocol (MCP) server for automated retrieval of scientific literature from PubMed and PubMed Central (PMC). The system provides programmatic access to search PubMed, identify relevant papers, and download full-text content or abstracts through a standardized interface compatible with AI assistant platforms.

### System Architecture

The PubMed Download MCP server was implemented as a Python-based service using the FastMCP framework (version ≥2.12.3), which provides Model Context Protocol compatibility for integration with AI systems. The architecture consists of three core modules:

1. **Server Module** (`server.py`): Implements the MCP protocol interface and exposes four tool endpoints: `search_pubmed_lr` (ligand-receptor search), `search_pubmed` (general search), `search_pubmed_minimal` (token-efficient search), and `download_paper` (paper retrieval).

2. **Search Module** (`utils/pubmed_search.py`): Handles all PubMed query construction and metadata retrieval using the NCBI Entrez API.

3. **Download Module** (`utils/paper_download.py`): Manages paper acquisition through multiple sources with intelligent fallback strategies.

All network operations were wrapped in asynchronous contexts using Python's `asyncio.to_thread()` to ensure non-blocking execution during MCP operations.

### Data Sources and APIs

The system interfaces with four primary data sources:

1. **NCBI Entrez E-utilities API**: Used for PubMed searches and metadata retrieval (accessed via Biopython ≥1.85).

2. **NCBI ID Converter API** (https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/): Used to check PMC availability for PubMed IDs.

3. **Jina Reader API** (https://r.jina.ai/): Used for full-text extraction from DOI and PMC URLs, with automatic conversion to clean markdown format.

4. **DOI Resolution System**: Used to resolve DOI identifiers to publisher URLs before content extraction.

### PubMed Search Implementation

#### General Search Method

The general search function accepts standard PubMed query syntax including field tags (e.g., `[Title]`, `[Author]`, `[Gene]`, `[PDAT]`), Boolean operators (`AND`, `OR`, `NOT`), and date filters. Search queries are submitted to the NCBI Entrez API using the `esearch` endpoint with the following parameters:

- Database: `pubmed`
- Return maximum: User-specified (default: 10 results)
- Sort order: `pub_date` (newest first) or `relevance`

After retrieving PubMed IDs (PMIDs), the system performs a batch fetch using `efetch` with XML return mode to retrieve complete metadata including title, abstract, authors, journal information, publication year, DOI, and PMC ID when available.

#### Ligand-Receptor Interaction Search

For ligand-receptor pair searches, queries are automatically constructed using the following template:

```
({ligand}[Gene Symbol] OR {ligand}[MeSH Terms]) AND
({receptor}[Gene Symbol] OR {receptor}[MeSH Terms]) AND
(interaction OR signaling OR binding OR pathway)
```

This query structure searches both Gene Symbol and MeSH Terms fields to maximize recall while focusing on papers discussing molecular interactions between the specified ligand and receptor.

#### Minimal Metadata Search

To optimize token usage for large-scale searches, a minimal search mode returns only essential fields: PMID, title, first author, publication year, journal, DOI, and PMC ID, omitting abstracts. This reduces token consumption by approximately 70-80% compared to full metadata retrieval.

### Paper Retrieval Strategy

#### Multi-Stage Fallback System

The paper download function implements a four-stage fallback strategy:

**Stage 1: Cache Check**
The system first checks for previously downloaded papers in three locations:
- Full-text from DOI: `./papers/{doi_sanitized}.txt`
- Full-text from PMID: `./papers/pmid_{pmid}.txt`
- Abstract-only: `./papers/pmid_{pmid}_abstract.txt`

If a cached version exists, it is returned immediately without network requests.

**Stage 2: DOI-Based Download**
When a DOI is available, the system:
1. Resolves the DOI to the publisher's URL using the https://doi.org/ redirect service
2. Constructs a Jina Reader API request: `https://r.jina.ai/{resolved_url}`
3. Submits the request with NCBI-compliant headers (see Compliance section)
4. Processes the response to extract clean text content

**Stage 3: PMC-Based Download**
When `prefer_pmc=True` and a PMID is provided:
1. Queries the NCBI ID Converter API to check PMC availability
2. If PMC ID is available, constructs PMC URL: `https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/`
3. Submits to Jina Reader API for content extraction
4. Falls back to DOI if PMC retrieval fails

**Stage 4: Abstract Fallback**
For papers without DOI or PMC availability (common for older publications):
1. Fetches metadata using Entrez `efetch` endpoint
2. Extracts title, authors, journal, year, PMID, and abstract
3. Formats as structured text file with metadata header
4. Saves to `./papers/pmid_{pmid}_abstract.txt`

Papers without abstracts return `[No abstract available]` rather than failing.

#### Download Retry Logic

All network requests implement exponential backoff retry with the following parameters:
- Maximum retries: 5 attempts
- Initial delay: 1.0 seconds
- Backoff multiplier: 2.0 (delays: 1s → 2s → 4s → 8s → 16s)
- Rate limiting detection: Automatically respects `Retry-After` headers for HTTP 429/503 responses

For persistent access failures (HTTP 401/403/405), the system progressively simplifies request headers and attempts search mode fallback (`https://s.jina.ai/`) to bypass cookie requirements or paywall restrictions.

### NCBI Compliance and Rate Limiting

To comply with NCBI usage policies and prevent server blocking:

1. **Email Identification**: All Entrez API requests include a user-provided email via `Entrez.email` setting, as required by NCBI.

2. **User-Agent Headers**: All requests use the identifier `PubMedDownloadMCP/1.0 ({email})` to clearly identify the client.

3. **Conservative Rate Limiting**:
   - Entrez API calls: 0.5 second delay between requests (2 req/s, below the 3 req/s limit)
   - PMC availability checks: 2-4 second random delay to prevent DDoS detection
   - DOI resolution: 1-2 second random delay
   - Optional NCBI API key support enables 10 req/s rate limit

4. **Request Spacing**: Random jitter is added to delays to avoid synchronized request patterns that could trigger rate limiting.

### Error Handling

The system implements comprehensive error handling:

- **Network Errors**: Wrapped in try-catch blocks with descriptive error messages returned to caller
- **API Failures**: Return structured dictionaries with `success: false` and detailed error descriptions
- **Timeout Protection**: All HTTP requests include 120-second timeouts to prevent hanging
- **Graceful Degradation**: PMC check failures silently fall back to DOI without raising exceptions
- **Cache Validation**: File existence checks before attempting downloads

### Data Output Formats

#### Search Results

All search functions return lists of dictionaries containing:
- `pmid`: PubMed ID (string)
- `title`: Article title (string)
- `abstract`: Abstract text (string, empty if unavailable)
- `authors`: Comma-separated author list (string)
- `year`: Publication year (string)
- `journal`: Journal title (string)
- `doi`: Digital Object Identifier (string, empty if unavailable)
- `pmcid`: PubMed Central ID (string, empty if unavailable)

Ligand-receptor searches additionally include `ligand` and `receptor` fields.

#### Download Results

The download function returns a dictionary with:
- `success`: Boolean indicating operation success
- `file_path`: Absolute path to downloaded file (string)
- `source`: Origin of content (`"Cache"`, `"DOI"`, `"PMC"`, or `"PubMed_Abstract"`)
- `message`: Human-readable status message (string)

#### File Formats

Full-text papers are saved as clean markdown-formatted text extracted by Jina Reader. Abstract-only files follow this format:

```
Title: {article_title}

Authors: {author_list}

Journal: {journal_name}
Year: {publication_year}
PMID: {pubmed_id}

Abstract:
{abstract_text}

---
Note: Full text not available. This file contains only the abstract and metadata from PubMed.
```

### Software Dependencies

The system requires Python ≥3.10 and the following packages:
- `fastmcp` ≥2.12.3 (MCP server framework)
- `biopython` ≥1.85 (NCBI Entrez API access)
- `requests` ≥2.32.5 (HTTP client)
- `bs4` ≥0.0.2 (HTML parsing for legacy support)

Development and debugging dependencies include:
- `ipykernel` ≥7.0.1 (interactive development)
- `pdbpp` ≥0.11.7 (enhanced debugger)
- `pykernel` ≥0.1.6 (kernel support)

### Configuration and Deployment

The server operates in standalone mode via stdio transport and requires environment variable configuration:
- `PUBMED_EMAIL` (required): Email for NCBI Entrez API identification
- `JINA_API_KEY` (optional): API key for enhanced Jina Reader access
- `NCBI_API_KEY` (optional): NCBI API key for increased rate limits (10 vs 3 req/s)

The server is launched via `uv run server.py` or `python server.py` and communicates through standard input/output using the Model Context Protocol specification.

### Reproducibility

All source code, documentation, and configuration examples are available in the project repository. The system maintains a flat package structure with core logic in `server.py`, `utils/pubmed_search.py`, and `utils/paper_download.py`. Complete implementation details, including rate limiting parameters, API endpoints, and error handling logic, are documented in the `CLAUDE.md` technical reference file included in the repository.
