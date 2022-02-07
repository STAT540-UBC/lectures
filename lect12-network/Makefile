RMD := $(wildcard *.rmd) $(wildcard *.Rmd)
PDF := $(RMD:.rmd=-slides.pdf) $(RMD:.Rmd=-slides.pdf)
HTML := $(RMD:.rmd=.html) $(RMD:.Rmd=.html)
TEX := $(RMD:.rmd=.tex) $(RMD:.Rmd=.tex)

all: $(HTML) $(PDF)

clean:
	rm -f $(PDF) $(HTML) $(TEX)

%.html: %.rmd
	Rscript -e "rmarkdown::render('$<', output_format='html_document', output_file='$@');"

%-slides.pdf: %.rmd %.html
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"

%.html: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='html_document', output_file='$@');"

%-slides.pdf: %.Rmd %.html
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"
