RMD := $(wildcard *.Rmd)
PDF := $(RMD:.Rmd=-slides.pdf)

all: $(PDF)

clean:
	rm -f $(PDF)

%-slides.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"
