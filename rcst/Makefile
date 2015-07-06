all: rcst.pdf

rcst.pdf: locate.pdf synth.pdf comp.pdf rcst.bbl rcst.tex
	pdflatex rcst
	pdflatex rcst

bib:
	pdflatex rcst
	bibtex rcst

rcst.bbl: rcst.bib
	pdflatex rcst
	bibtex rcst

locate.pdf: locate.csv locate.R
	R --slave --args locate xy legend < locate.R

synth.pdf: synth.csv synth.R
	R --slave --args synth xy legend < synth.R

comp.pdf: comp.csv synth.R
	R --slave --args comp x legend < synth.R

clean:
	rm -f locate.pdf synth.pdf comp.pdf rcst.aux rcst.bbl rcst.bib.bak rcst.blg rcst.log rcst.pdf
