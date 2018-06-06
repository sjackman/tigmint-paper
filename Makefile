pandoc_opt=-Fpandoc-crossref -Fpandoc-citeproc

.DELETE_ON_ERROR:
.SECONDARY:

all: tigmint.html tigmint.pdf

docx: tigmint.docx

clean:
	rm -f tigmint.html tigmint.pdf tigmint-supp.html tigmint-supp.pdf

# Download the citation style language (CSL).
tigmint.csl:
	curl -o $@ https://www.zotero.org/styles/biomed-central

# Render Markdown to HTML using Pandoc.
%.html: %.md
	pandoc $(pandoc_opt) -s --mathjax -o $@ $<

# Render Markdown to PDF using Pandoc.
%.pdf: %.md
	pandoc $(pandoc_opt) -o $@ $<

# Render Markdown to DOCX using Pandoc.
%.docx: %.md
	pandoc -o $@ $<

# Generate Table of Contents for supplemental material only
tigmint-supp.pdf: tigmint-supp.md
	pandoc $(pandoc_opt) --toc -o $@ $<

# Fetch BibTex records from a list of DOI.
%.doi.bib: %.doi
	brew cite $$(<$<) | sort >$@

# Concatentate the citations with and without DOI.
%.bib: %.doi.bib %.nodoi.bib
	sort $^ | sed 's~http://dx.doi.org~https://doi.org~' >$@

tigmint.docx: tigmint.bib tigmint.csl
tigmint.html: tigmint.bib tigmint.csl
tigmint.pdf: tigmint.bib tigmint.csl

tigmint-supp.docx: tigmint.bib tigmint.csl
tigmint-supp.html: tigmint.bib tigmint.csl
tigmint-supp.pdf: tigmint.bib tigmint.csl
