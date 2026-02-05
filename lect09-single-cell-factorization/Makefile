## Find .Rmd files (beamer slides) and .rmd files (regular PDF documents)
RMD_SLIDES := $(wildcard *.Rmd)

## Child documents and dependencies
CHILD_RMD := $(wildcard lect09-*-section.rmd) $(wildcard lect09-*-applications.rmd)
TEXFRAG := $(wildcard texfrag/*.tex)

## Output targets
PDF_SLIDES := $(RMD_SLIDES:.Rmd=-slides.pdf)

# HTML := $(RMD_DOCS:.rmd=.html) $(RMD_SLIDES:.Rmd=.html)
# TEX := $(RMD_DOCS:.rmd=.tex) $(RMD_SLIDES:.Rmd=.tex)

all: $(PDF_SLIDES) $(PDF_DOCS)

clean:
	rm -f $(PDF_SLIDES) $(PDF_DOCS) $(RMD_SLIDES:.Rmd=-handout.pdf)
	rm -f *.log *.tex

## Specific rule for lect09 with dependencies on child documents and texfrag
lect09-single-cell-factorization-slides.pdf: lect09-single-cell-factorization.Rmd $(CHILD_RMD) $(TEXFRAG) Setup.R
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"
	rm -f lect09-single-cell-factorization-slides.log lect09-single-cell-factorization-slides.tex

## Beamer slides from .Rmd files (generic rule)
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
