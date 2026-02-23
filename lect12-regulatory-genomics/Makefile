## Find .Rmd files (beamer slides) and .rmd files (regular PDF documents)
RMD_SLIDES := $(wildcard *.Rmd)
RMD_DOCS := $(wildcard *.rmd)

## Output targets
PDF_SLIDES := $(RMD_SLIDES:.Rmd=-slides.pdf)
PDF_DOCS := $(RMD_DOCS:.rmd=.pdf)

HTML := $(RMD_DOCS:.rmd=.html) $(RMD_SLIDES:.Rmd=.html)
TEX := $(RMD_DOCS:.rmd=.tex) $(RMD_SLIDES:.Rmd=.tex)

all: $(PDF_SLIDES) $(PDF_DOCS)

clean:
	rm -f $(PDF_SLIDES) $(PDF_DOCS)
	rm -f *.log *.tex

## Beamer slides from .Rmd files
%-slides.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"
	rm -f $*-slides.log $*-slides.tex

## Handout from .Rmd files
%-handout.pdf: %.Rmd handout-header.tex
	Rscript -e "rmarkdown::render('$<', output_format=rmarkdown::beamer_presentation(pandoc_args=c('--variable=classoption:handout,aspectratio=169', '--include-in-header=handout-header.tex')), output_file='$@');"
	rm -f $*-handout.log $*-handout.tex

## Regular PDF documents from .rmd files
%.pdf: %.rmd
	Rscript -e "rmarkdown::render('$<', output_format='pdf_document', output_file='$@');"
	rm -f $*.log $*.tex
