all:
	#"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit2html('ChemmineR.Rmd',options=c('toc',markdown::markdownHTMLOptions(TRUE)))"
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit('ChemmineR.Rmd')"
	pandoc --output ChemmineR.html -H header.html --from markdown --smart --table-of-contents ChemmineR.md;

clean:
	rm -f *.md
