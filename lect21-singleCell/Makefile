RMD := $(wildcard *.rmd) $(wildcard *.Rmd)
PDF := $(RMD:.rmd=-slides.pdf) $(RMD:.Rmd=-slides.pdf)
HTML := $(RMD:.rmd=.html) $(RMD:.Rmd=.html)
TEX := $(RMD:.rmd=.tex) $(RMD:.Rmd=.tex)

all: $(PDF)

clean:
	rm -f $(PDF)

%-slides.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"
