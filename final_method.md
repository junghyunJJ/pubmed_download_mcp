# Paper Retrieval Strategy - Final Method

To support vector database construction and Retrieval-Augmented Generation (RAG) applications, we developed a paper retrieval system that implements a four-stage fallback strategy to maximize retrieval success while minimizing redundant network requests. First, the system checks for previously downloaded papers in local cache locations, including full-text files and abstract-only files, returning cached versions immediately when available. Second, when a DOI is provided, the system resolves it to the publisher's URL using the DOI resolution system [8], constructs a Jina Reader API [7] request, and extracts clean text content. Third, for papers with PMIDs when `prefer_pmc=True` is specified, the system queries the NCBI ID Converter API [6] to determine PMC availability and, if present, retrieves full text from PubMed Central, falling back to DOI-based retrieval if PMC access fails. Fourth, for older publications lacking both DOI and PMC availability, the system fetches title, authors, journal, year, and abstract metadata using the Entrez `efetch` endpoint [4] and saves the structured information to abstract-only files, returning "[No abstract available]" for papers without abstracts rather than failing. Additionally, paper searches were conducted using bioMCP [13], which provides access to PubMed and PubTator3 databases for comprehensive literature search and biomedical entity annotation. Once relevant articles were identified through bioMCP, the full-text papers or abstracts were downloaded using the retrieval strategy described above.

## References

[4] National Center for Biotechnology Information. Entrez Programming Utilities Help. Bethesda (MD): National Library of Medicine (US); 2010. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25501/

[6] National Center for Biotechnology Information. PMC ID Converter API. Available from: https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/

[7] Jina AI. Jina Reader API Documentation. Available from: https://jina.ai/reader/

[8] International DOI Foundation. DOI Resolution System. Available from: https://www.doi.org/

[13] bioMCP. Finding Articles and cBioPortal Data. Available from: https://biomcp.org/how-to-guides/01-find-articles-and-cbioportal-data/
