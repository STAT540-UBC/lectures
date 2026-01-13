RMD := $(wildcard *.Rmd)
SLIDES := $(RMD:.Rmd=-slides.pdf)
HANDOUT := $(RMD:.Rmd=-handout.pdf)
HTML := $(RMD:.rmd=.html) $(RMD:.Rmd=.html)
TEX := $(RMD:.rmd=.tex) $(RMD:.Rmd=.tex)

all: $(SLIDES)

handout: $(HANDOUT)

clean:
	rm -f $(SLIDES) $(HANDOUT)

%-slides.pdf: %.Rmd
	Rscript -e "rmarkdown::render('$<', output_format='beamer_presentation', output_file='$@');"
	rm -f $*-slides.log

%-handout.pdf: %.Rmd handout-header.tex
	Rscript -e "rmarkdown::render('$<', output_format=rmarkdown::beamer_presentation(pandoc_args=c('--variable=classoption:handout,aspectratio=169', '--include-in-header=handout-header.tex')), output_file='$@');"
	rm -f %-handout.log  %-handout.tex

