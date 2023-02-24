RMD := $(wildcard *.Rmd)
PDF := $(RMD:.Rmd=-slides.pdf) $(RMD:.Rmd=-slides.pdf)
HTML := $(RMD:.Rmd=.html) $(RMD:.Rmd=.html)
TEX := $(RMD:.Rmd=.tex) $(RMD:.Rmd=.tex)

all: $(PDF)

clean:
	rm -f $(PDF)

%-slides.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format=rmarkdown::beamer_presentation(fig_crop=FALSE), output_file='$@');"
