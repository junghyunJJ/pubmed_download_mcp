# Paper Retrieval Strategy

To support vector database construction and Retrieval-Augmented Generation (RAG) applications [12], we developed a paper retrieval system that implements a four-stage fallback strategy to maximize retrieval success while minimizing redundant network requests. First, the system checks for previously downloaded papers in local cache locations, including full-text files (`./papers/{doi_sanitized}.txt` and `./papers/pmid_{pmid}.txt`) and abstract-only files (`./papers/pmid_{pmid}_abstract.txt`), returning cached versions immediately when available. Second, when a DOI is provided, the system resolves it to the publisher's URL using the DOI resolution system [8], constructs a Jina Reader API request [7], and extracts clean text content. Third, for papers with PMIDs when `prefer_pmc=True` is specified, the system queries the NCBI ID Converter API [6] to determine PMC availability and, if present, retrieves full text from PubMed Central, falling back to DOI-based retrieval if PMC access fails. Fourth, for older publications lacking both DOI and PMC availability, the system fetches title, authors, journal, year, and abstract metadata using the Entrez `efetch` endpoint [4] and saves the structured information to abstract-only files, returning "[No abstract available]" for papers without abstracts rather than failing.

All network requests implement exponential backoff retry with up to 5 attempts, initial delay of 1.0 seconds, and backoff multiplier of 2.0 (delays: 1s → 2s → 4s → 8s → 16s) [9]. The system automatically respects `Retry-After` headers for HTTP 429/503 responses and progressively simplifies request headers for persistent access failures (HTTP 401/403/405), ultimately attempting search mode fallback via `https://s.jina.ai/` [7] to bypass authentication requirements. To comply with NCBI usage policies [10] and prevent rate limiting, the system implements conservative delays including 0.5 seconds between Entrez API calls (2 requests/second, below the 3 requests/second limit), 2-4 second random delays for PMC availability checks, and 1-2 second random delays for DOI resolution. All requests include NCBI-compliant User-Agent headers (`PubMedDownloadMCP/1.0 {email}`) with user-provided email addresses [4,10], and optional NCBI API key support enables increased rate limits of 10 requests/second [11]. Error handling is implemented through try-catch blocks with descriptive messages, 120-second HTTP timeouts to prevent hanging [9], and graceful degradation where PMC check failures silently fall back to alternative retrieval methods without raising exceptions.

Paper searches were conducted using bioMCP's Finding Articles feature [13], which provides access to PubMed and PubTator3 databases for comprehensive literature search and biomedical entity annotation. Once relevant articles were identified through bioMCP, the full-text papers or abstracts were downloaded using the retrieval strategy described above.

## References

[4] National Center for Biotechnology Information. Entrez Programming Utilities Help. Bethesda (MD): National Library of Medicine (US); 2010. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25501/

[6] National Center for Biotechnology Information. PMC ID Converter API. Available from: https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

[7] Jina AI. Jina Reader API Documentation. Available from: https://jina.ai/reader/

[8] International DOI Foundation. DOI Resolution System. Available from: https://www.doi.org/

[9] Reitz K. Requests: HTTP for Humans. Python Software Foundation; 2020. Available from: https://requests.readthedocs.io/

[10] National Center for Biotechnology Information. E-utilities Usage Guidelines and Requirements. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen

[11] NCBI Insights. New API Keys for the E-utilities. 2017 Nov 2. Available from: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/

[12] Lewis P, Perez E, Piktus A, Petroni F, Karpukhin V, Goyal N, et al. Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks. In: Advances in Neural Information Processing Systems 33 (NeurIPS 2020); 2020. p. 9459-9474.

[13] bioMCP. Finding Articles and cBioPortal Data. Available from: https://biomcp.org/how-to-guides/01-find-articles-and-cbioportal-data/
