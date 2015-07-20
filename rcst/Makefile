SOURCES=$(wildcard *.tex)

all: paper.pdf

paper.pdf: locate.pdf synth.pdf comp.pdf paper.bbl $(SOURCES)
	pdflatex paper
	pdflatex paper

bib:
	pdflatex paper
	bibtex paper

paper.bbl: paper.bib
	pdflatex paper
	bibtex paper

locate.pdf: locate.csv locate.R
	R --slave --args locate xy legend < locate.R

synth.pdf: synth.csv synth.R
	R --slave --args synth xy legend < synth.R

comp.pdf: comp.csv synth.R
	R --slave --args comp x legend < synth.R

clean:
	rm -f *.aux *.bak *.bbl *.blg *.log *.pdf
