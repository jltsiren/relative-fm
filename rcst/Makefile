SOURCES=$(wildcard *.tex)
FIGURES=locate.pdf locate_bars.pdf synth.pdf comp.pdf

all: paper.pdf

paper.pdf: $(SOURCES) $(FIGURES) paper.bbl
	pdflatex paper
	pdflatex paper

force: $(FIGURES) paper.bbl
	rm -f paper.pdf
	pdflatex paper
	pdflatex paper

bib: $(FIGURES)
	pdflatex paper
	bibtex paper

paper.bbl: paper.bib
	pdflatex paper
	bibtex paper

locate.pdf: locate.csv locate.R
	R --slave --args locate x legend < locate.R

locate_bars.pdf: locate.csv locate_bars.R
	R --slave --args locate xy legend < locate_bars.R

synth.pdf: synth.csv synth.R
	R --slave --args synth xy legend < synth.R

comp.pdf: comp.csv synth.R
	R --slave --args comp x legend < synth.R

clean:
	rm -f *.aux *.bak *.bbl *.blg *.log paper.pdf $(FIGURES)
