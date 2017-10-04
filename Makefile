pandoc_opt=-Fpandoc-citeproc -Fpandoc-crossref

.DELETE_ON_ERROR:
.SECONDARY:

all: tigmint.html

clean:
	rm -f tigmint.html tigmint.pdf tigmint-supp.html tigmint-supp.pdf

# Download the citation style language (CSL).
tigmint.csl:
	curl -o $@ https://www.zotero.org/styles/bioinformatics

# Render Markdown to HTML using Pandoc.
%.html: %.md
	pandoc $(pandoc_opt) -s -o $@ $<

# Render Markdown to PDF using Pandoc.
%.pdf: %.md
	pandoc $(pandoc_opt) -o $@ $<

# Generate Table of Contents for supplemental material only
tigmint-supp.pdf: tigmint-supp.md
	pandoc $(pandoc_opt) --toc -o $@ $<

# Render Markdown to DOCX using Pandoc.
%.docx: %.md
	pandoc $(pandoc_opt) -o $@ $<

# Render RMarkdown to HTML using R.
%.html: %.rmd
	RScript -e 'rmarkdown::render("$<")'

tigmint.docx: tigmint.bib tigmint.csl
tigmint.html: tigmint.bib tigmint.csl
tigmint.pdf: tigmint.bib tigmint.csl

tigmint-supp.docx: tigmint.bib tigmint.csl
tigmint-supp.html: tigmint.bib tigmint.csl
tigmint-supp.pdf: tigmint.bib tigmint.csl
