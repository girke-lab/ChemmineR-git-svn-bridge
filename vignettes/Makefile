all:
	"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e "library(knitr); knit2html('ChemmineR.Rmd',options=c('toc',markdown::markdownHTMLOptions(TRUE)))"
#	pandoc -s ChemmineR.html -o ChemmineR.pdf
clean:
	rm -f *.md
