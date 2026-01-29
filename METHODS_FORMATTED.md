# Methods

## Development of PubMed Literature Retrieval System

### Overview

We developed a Model Context Protocol (MCP) server for automated retrieval of scientific literature from PubMed and PubMed Central (PMC). The system provides programmatic access to search PubMed, identify relevant papers, and download full-text content or abstracts through a standardized interface compatible with AI assistant platforms. The server was implemented as a Python-based service [1] using the FastMCP framework (version ≥2.12.3) [2], which provides Model Context Protocol compatibility [3] for integration with AI systems.

### System Architecture

The PubMed Download MCP server consists of three core modules. The server module (`server.py`) implements the MCP protocol interface [3] and exposes four tool endpoints: `search_pubmed_lr` (ligand-receptor search), `search_pubmed` (general search), `search_pubmed_minimal` (token-efficient search), and `download_paper` (paper retrieval). The search module (`utils/pubmed_search.py`) handles all PubMed query construction and metadata retrieval using the NCBI Entrez API [4]. The download module (`utils/paper_download.py`) manages paper acquisition through multiple sources with intelligent fallback strategies. All network operations were wrapped in asynchronous contexts using Python's `asyncio.to_thread()` [1] to ensure non-blocking execution during MCP operations.

### Data Sources and APIs

The system interfaces with four primary data sources. PubMed searches and metadata retrieval were performed using the NCBI Entrez E-utilities API [4], accessed via Biopython (version ≥1.85) [5]. PMC availability checks utilized the NCBI ID Converter API [6]. Full-text extraction from DOI and PMC URLs employed the Jina Reader API [7], which provides automatic conversion to clean markdown format. DOI identifiers were resolved to publisher URLs using the DOI resolution system [8] before content extraction.

### PubMed Search Implementation

#### General Search Method

The general search function accepts standard PubMed query syntax including field tags (e.g., `[Title]`, `[Author]`, `[Gene]`, `[PDAT]`), Boolean operators (`AND`, `OR`, `NOT`), and date filters [4]. Search queries were submitted to the NCBI Entrez API using the `esearch` endpoint with database parameter set to `pubmed`, user-specified return maximum (default: 10 results), and sort order by publication date (newest first) or relevance. After retrieving PubMed IDs (PMIDs), the system performed batch fetches using `efetch` with XML return mode to retrieve complete metadata including title, abstract, authors, journal information, publication year, DOI, and PMC ID when available [4].

#### Ligand-Receptor Interaction Search

For ligand-receptor pair searches, queries were automatically constructed using the following template: `({ligand}[Gene Symbol] OR {ligand}[MeSH Terms]) AND ({receptor}[Gene Symbol] OR {receptor}[MeSH Terms]) AND (interaction OR signaling OR binding OR pathway)`. This query structure searches both Gene Symbol and MeSH Terms fields to maximize recall while focusing on papers discussing molecular interactions between the specified ligand and receptor.

#### Minimal Metadata Search

To optimize token usage for large-scale searches, a minimal search mode returns only essential fields: PMID, title, first author, publication year, journal, DOI, and PMC ID, omitting abstracts. This approach reduces token consumption by approximately 70-80% compared to full metadata retrieval.

### Paper Retrieval Strategy

#### Multi-Stage Fallback System

The paper download function implements a four-stage fallback strategy. First, the system checks for previously downloaded papers in three cache locations: full-text from DOI (`./papers/{doi_sanitized}.txt`), full-text from PMID (`./papers/pmid_{pmid}.txt`), and abstract-only (`./papers/pmid_{pmid}_abstract.txt`). If a cached version exists, it is returned immediately without network requests.

When a DOI is available, the system resolves the DOI to the publisher's URL using the https://doi.org/ redirect service [8], constructs a Jina Reader API request (`https://r.jina.ai/{resolved_url}`) [7], and processes the response to extract clean text content. For papers with PMIDs when `prefer_pmc=True`, the system queries the NCBI ID Converter API [6] to check PMC availability and, if available, submits the PMC URL to the Jina Reader API for content extraction, falling back to DOI if PMC retrieval fails.

For papers without DOI or PMC availability (common for older publications), the system fetches metadata using the Entrez `efetch` endpoint [4], extracts title, authors, journal, year, PMID, and abstract, formats the information as structured text with metadata header, and saves to `./papers/pmid_{pmid}_abstract.txt`. Papers without abstracts return "[No abstract available]" rather than failing.

#### Download Retry Logic

All network requests implement exponential backoff retry with maximum 5 attempts, initial delay of 1.0 seconds, backoff multiplier of 2.0 (delays: 1s → 2s → 4s → 8s → 16s), and automatic respect for `Retry-After` headers for HTTP 429/503 responses [9]. For persistent access failures (HTTP 401/403/405), the system progressively simplifies request headers and attempts search mode fallback (`https://s.jina.ai/`) [7] to bypass cookie requirements or paywall restrictions.

### NCBI Compliance and Rate Limiting

To comply with NCBI usage policies [10] and prevent server blocking, all Entrez API requests include a user-provided email via `Entrez.email` setting [4,5]. All requests use the User-Agent identifier `PubMedDownloadMCP/1.0 ({email})` to clearly identify the client. Conservative rate limiting includes 0.5 second delays between Entrez API calls (2 requests/second, below the 3 requests/second limit) [10], 2-4 second random delays for PMC availability checks to prevent DDoS detection, 1-2 second random delays for DOI resolution, and optional NCBI API key support enabling 10 requests/second rate limit [11]. Random jitter is added to delays to avoid synchronized request patterns that could trigger rate limiting.

### Error Handling

The system implements comprehensive error handling with network errors wrapped in try-catch blocks with descriptive error messages returned to the caller, API failures returning structured dictionaries with `success: false` and detailed error descriptions, 120-second timeouts for all HTTP requests [9] to prevent hanging, graceful degradation where PMC check failures silently fall back to DOI without raising exceptions, and cache validation with file existence checks before attempting downloads.

### Data Output Formats

#### Search Results

All search functions return lists of dictionaries containing: `pmid` (PubMed ID as string), `title` (article title as string), `abstract` (abstract text as string, empty if unavailable), `authors` (comma-separated author list as string), `year` (publication year as string), `journal` (journal title as string), `doi` (Digital Object Identifier as string, empty if unavailable), and `pmcid` (PubMed Central ID as string, empty if unavailable). Ligand-receptor searches additionally include `ligand` and `receptor` fields.

#### Download Results

The download function returns a dictionary with: `success` (Boolean indicating operation success), `file_path` (absolute path to downloaded file as string), `source` (origin of content: `"Cache"`, `"DOI"`, `"PMC"`, or `"PubMed_Abstract"`), and `message` (human-readable status message as string).

#### File Formats

Full-text papers are saved as clean markdown-formatted text extracted by Jina Reader [7]. Abstract-only files follow a structured format with title, authors, journal, year, PMID, abstract text, and a note indicating that full text was not available.

### Software Dependencies

The system requires Python ≥3.10 [1] and the following packages: `fastmcp` ≥2.12.3 (MCP server framework) [2], `biopython` ≥1.85 (NCBI Entrez API access) [5], `requests` ≥2.32.5 (HTTP client) [9], and `bs4` ≥0.0.2 (HTML parsing for legacy support) [12]. Development and debugging dependencies include `ipykernel` ≥7.0.1 (interactive development), `pdbpp` ≥0.11.7 (enhanced debugger), and `pykernel` ≥0.1.6 (kernel support).

### Configuration and Deployment

The server operates in standalone mode via stdio transport and requires environment variable configuration: `PUBMED_EMAIL` (required, email for NCBI Entrez API identification) [4,10], `JINA_API_KEY` (optional, API key for enhanced Jina Reader access) [7], and `NCBI_API_KEY` (optional, NCBI API key for increased rate limits from 3 to 10 requests/second) [11]. The server is launched via `uv run server.py` or `python server.py` and communicates through standard input/output using the Model Context Protocol specification [3].

### Reproducibility

All source code, documentation, and configuration examples are available in the project repository. The system maintains a flat package structure with core logic in `server.py`, `utils/pubmed_search.py`, and `utils/paper_download.py`. Complete implementation details, including rate limiting parameters, API endpoints, and error handling logic, are documented in the technical reference file included in the repository.

## References

[1] Van Rossum G, Drake FL. Python 3 Reference Manual. Scotts Valley, CA: CreateSpace; 2009.

[2] FastMCP: Model Context Protocol Server Framework. Available from: https://github.com/jlowin/fastmcp

[3] Anthropic. Model Context Protocol Specification. Available from: https://modelcontextprotocol.io/

[4] National Center for Biotechnology Information. Entrez Programming Utilities Help. Bethesda (MD): National Library of Medicine (US); 2010. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25501/

[5] Cock PJ, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25(11):1422-1423.

[6] National Center for Biotechnology Information. PMC ID Converter API. Available from: https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

[7] Jina AI. Jina Reader API Documentation. Available from: https://jina.ai/reader/

[8] International DOI Foundation. DOI Resolution System. Available from: https://www.doi.org/

[9] Reitz K. Requests: HTTP for Humans. Python Software Foundation; 2020. Available from: https://requests.readthedocs.io/

[10] National Center for Biotechnology Information. E-utilities Usage Guidelines and Requirements. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen

[11] NCBI Insights. New API Keys for the E-utilities. 2017 Nov 2. Available from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

[12] Richardson L. Beautiful Soup Documentation. Available from: https://www.crummy.com/software/BeautifulSoup/bs4/doc/
