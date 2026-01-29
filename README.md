# PubMed Download MCP Server

A Model Context Protocol (MCP) server for searching PubMed and downloading scientific papers. Built with FastMCP for seamless integration with Claude Desktop.

## Features

- **PubMed Search**: Search papers by ligand-receptor pairs or general queries
- **Paper Download**: Download full-text papers from PMC (via NCBI efetch) or DOI
- **NCBI efetch**: Official NCBI API for PMC downloads (no DDoS blocking)
- **JATS XML Parsing**: Clean plain text extraction from PMC's XML format
- **Abstract Fallback**: Automatically fetch abstracts for papers without full text
- **Smart Caching**: Avoid redundant downloads with automatic caching
- **Async Operations**: Non-blocking operations for better performance
- **NCBI Compliant**: Rate limiting and proper headers to prevent blocking

## Version 4.1 Updates

### New Features (v4.1)
- **Improved Fallback Chain**: PMC failure now properly falls back to DOI
- **Modular Code**: Helper functions for cleaner, maintainable code
- **Better Error Handling**: Detailed failure messages showing all attempted methods

### Previous Features (v4.0)
- **NCBI efetch for PMC**: Direct download via official NCBI API (no DDoS blocking)
- **JATS XML Parsing**: Full-text extraction from PMC's XML format
- **References Excluded**: Cleaner output for RAG (removes reference noise)
- **PubMed Abstract Fallback**: Papers without PMC or DOI return abstracts with metadata

### Download Strategy (v4.1)
```
1. Check Cache (full text or abstract)
2. If PMC available: Try NCBI efetch (official API, no blocking)
3. If efetch fails: Try PMC via Jina Reader
4. If PMC fails: Try DOI via Jina Reader  ← NEW in v4.1
5. Final Fallback: Fetch abstract from PubMed
```

## Installation

### Prerequisites
- Python 3.8+
- uv (recommended) or pip
- NCBI Entrez email (required for PubMed API)
- Jina API key (optional, for faster downloads)

### Install Dependencies

```bash
# Using uv (recommended)
uv sync

# Or using pip
pip install -e .
```

## Configuration

### Environment Variables

Create a `.env` file or set environment variables:

```bash
# Required
export PUBMED_EMAIL="your-email@example.com"

# Optional
export JINA_API_KEY="your-jina-api-key"
export NCBI_API_KEY="your-ncbi-api-key"  # For 10 req/s rate limit vs 3 req/s
```

### Claude Desktop Configuration

Add to `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "pubmed-download": {
      "command": "uv",
      "args": [
        "--directory",
        "/absolute/path/to/pubmed_download_mcp",
        "run",
        "server.py"
      ],
      "env": {
        "PUBMED_EMAIL": "your-email@example.com",
        "JINA_API_KEY": "your-optional-api-key",
        "NCBI_API_KEY": "your-optional-ncbi-key"
      }
    }
  }
}
```

## Usage

### Available Tools

#### 1. `search_pubmed_lr`
Search for ligand-receptor interaction papers.

```python
search_pubmed_lr(
    ligand="IL2",
    receptor="IL2RA",
    max_results=5,
    sort_by="pub_date"  # or "relevance"
)
```

**Returns**: List of papers with PMID, title, abstract, authors, year, journal, DOI, PMC ID

#### 2. `search_pubmed`
General PubMed search with full query syntax support.

```python
search_pubmed(
    query="BRCA1[Gene] AND breast cancer",
    max_results=10,
    sort_by="pub_date"
)
```

**Supports**: Field tags (`[Title]`, `[Author]`, `[Gene]`), Boolean operators (`AND`, `OR`, `NOT`)

#### 3. `search_pubmed_minimal`
Token-efficient search returning only essential metadata (no abstracts).

```python
search_pubmed_minimal(
    query="cancer immunotherapy",
    max_results=10
)
```

**Returns**: PMID, title, first author, year, journal, DOI, PMC ID (70-80% fewer tokens)

#### 4. `download_paper`
Download papers with smart fallback strategy.

```python
download_paper(
    pmid="33239918",           # Optional if DOI provided
    doi="10.1038/...",         # Optional if PMID provided
    output_path=None,          # Auto-generates if not provided
    prefer_pmc=True            # Try PMC before DOI
)
```

**Returns**: Dictionary with success status, file path, source, and message

### Download Behavior (v4.1)

| Scenario | Source | Result |
|----------|--------|--------|
| Already cached | `Cache` | Cached file (no download) |
| PMC available + efetch works | `PMC_efetch` | Full-text (JATS XML → plain text) |
| PMC efetch fails + Jina works | `PMC` | Full-text (via Jina Reader) |
| PMC fails + DOI available | `DOI` | Full-text (via Jina Reader) ← **v4.1** |
| No PMC/DOI available | `PubMed_Abstract` | Abstract + metadata |

### Output Files

```
papers/
├── pmid_33239918.txt              # Full-text from PMC
├── 10.1038_nature12345.txt        # Full-text from DOI
└── pmid_10492251_abstract.txt     # Abstract only (no PMC/DOI)
```

## Examples

### Example 1: Search and Download Recent Papers

```python
# Search for recent papers
results = await search_pubmed(
    query="prostate cancer[Title] AND 2024[PDAT]",
    max_results=5
)

# Download first paper
for paper in results:
    result = await download_paper(
        pmid=paper["pmid"],
        prefer_pmc=True
    )
    print(f"Downloaded: {result['file_path']}")
```

### Example 2: Download with Fallback

```python
# Try PMID 10492251 (old paper, no PMC/DOI)
result = await download_paper(pmid="10492251")

# Result:
# {
#   "success": True,
#   "file_path": "./papers/pmid_10492251_abstract.txt",
#   "source": "PubMed_Abstract",
#   "message": "Successfully saved abstract to ..."
# }
```

### Example 3: Ligand-Receptor Search

```python
# Search for IL-2/IL-2RA interactions
papers = await search_pubmed_lr(
    ligand="IL2",
    receptor="IL2RA",
    max_results=10
)

# Returns papers discussing IL-2 and IL-2RA interaction, signaling, binding
```

## File Formats

### Full-Text Papers (PMC efetch)
- Clean plain text parsed from JATS XML
- Includes: title, authors, journal, year, DOI, PMC ID
- Includes: abstract (handles structured abstracts)
- Includes: full body text with section headers
- References excluded (for RAG efficiency)

### Full-Text Papers (Jina Reader)
- Clean markdown format from Jina Reader
- Includes title, authors, abstract, full body text

### Abstract-Only Files
```
Title: [Paper title]

Authors: [Author list]

Journal: [Journal name]
Year: [Year]
PMID: [PMID]

Abstract:
[Abstract text]

---
Note: Full text not available. This file contains only the abstract and metadata from PubMed.
```

## Rate Limiting

The server implements NCBI-compliant rate limiting:

- **PubMed Search**: 0.5s delay between requests (2 req/s)
- **PMC efetch**: 0.5s delay (official API, no DDoS issues)
- **PMC Checks**: 2-4s random delay (to prevent DDoS detection)
- **DOI Requests**: 1-2s delay
- **With NCBI API Key**: Up to 10 req/s (vs 3 req/s without)

## Troubleshooting

### Common Issues

**Issue**: "PUBMED_EMAIL environment variable required"
- **Solution**: Set `PUBMED_EMAIL` in environment or Claude Desktop config

**Issue**: "Paper not available"
- **Check**: Run with just PMID - will fall back to abstract if full text unavailable
- **Note**: Some old papers (pre-2000) may not have DOI or PMC availability

**Issue**: Papers downloading slowly
- **Solution**: Add `JINA_API_KEY` for faster Jina Reader access
- **Solution**: Add `NCBI_API_KEY` for higher rate limits (10 req/s)

**Issue**: Rate limiting errors
- **Cause**: Too many requests to NCBI/PMC
- **Solution**: Server automatically implements delays; wait and retry

**Issue**: Jina Reader blocked by PMC (DDoS detection)
- **Solution**: v4.0+ uses NCBI efetch as primary method (no blocking)
- **Note**: efetch only works for PMC Open Access papers
- **v4.1**: If PMC fails completely, now falls back to DOI automatically

### Debug Mode

Run server directly to see detailed logs:

```bash
export PUBMED_EMAIL="your-email@example.com"
uv run server.py
```

Check logs:
```bash
# Project log
tail -f pubmed_mcp.log

# Claude Desktop logs (macOS)
tail -f ~/Library/Logs/Claude/mcp-server-pubmed-download.log
```

## API Keys

### NCBI API Key
- **Get**: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
- **Benefit**: 10 requests/second vs 3 requests/second
- **Required**: No (optional)

### Jina API Key
- **Get**: https://jina.ai/
- **Benefit**: Faster paper downloads, better reliability
- **Required**: No (optional)

## Development

### Running Tests

```bash
# Test paper download
uv run test_download.py

# Test specific PMID
python -c "
import asyncio
from utils.paper_download import download_paper_from_doi

async def test():
    result = await download_paper_from_doi(pmid='33239918')
    print(result)

asyncio.run(test())
"
```

### Project Structure

```
pubmed_download_mcp/
├── server.py                    # Main MCP server
├── utils/
│   ├── pubmed_search.py        # PubMed search functions
│   ├── paper_download.py       # Paper download logic
│   └── pmc_to_lmm_txt.py      # Deprecated PMC HTML parser
├── papers/                      # Downloaded papers directory
├── pyproject.toml              # Dependencies
└── README.md                    # This file
```

## Contributing

Issues and pull requests are welcome! Please ensure:
- Code follows existing patterns (async wrappers, error handling)
- NCBI rate limiting is respected
- Changes maintain backward compatibility

## License

MIT License - See LICENSE file for details

## Changelog

### v4.1 (Current)
- ✨ **Improved fallback chain**: PMC failure now properly falls back to DOI
- ✨ Added `_try_jina_download()` helper for modular Jina API calls
- ✨ Added `_determine_output_path()` and `_save_content()` helpers
- ✨ Better error messages showing all attempted download methods
- 🐛 Fixed bug where PMC Jina failure skipped DOI fallback

### v4.0
- ✨ Added NCBI efetch for PMC downloads (official API, no DDoS blocking)
- ✨ Added JATS XML to plain text parser (`parse_jats_xml_to_text`)
- ✨ References excluded from output for RAG efficiency
- ✨ New source type: `PMC_efetch`
- 🔧 Added `lxml` dependency for XML parsing

### v3.0
- ✨ Added PubMed abstract fallback for papers without PMC/DOI
- ✨ Enhanced cache checking (full text + abstract)
- ✨ Removed `force_download` option from MCP tool
- ✨ Improved error messages and graceful degradation
- 🐛 Fixed handling of papers without abstracts

### v2.1
- 🔧 Simplified download strategy (removed PMC HTML parser)
- 🔧 Direct Jina Reader usage for PMC URLs
- 🐛 Fixed 403 errors with NCBI-compliant headers

### v2.0
- ✨ Added support for PMID-only downloads
- ✨ Added `prefer_pmc` option (default: True)
- 🔧 Improved PMC availability checking with NCBI ID Converter API

### v1.0
- 🎉 Initial release with basic search and download functionality

## Support

For issues or questions:
1. Check the troubleshooting section above
2. Review CLAUDE.md for detailed technical documentation
3. Check Claude Desktop logs for error details
4. File an issue with reproduction steps and error messages
