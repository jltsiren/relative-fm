all: rcst.pdf

rcst.pdf: locate.pdf rcst.bbl rcst.tex
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

clean:
	rm -f locate.pdf rcst.aux rcst.bbl rcst.bib.bak rcst.blg rcst.log rcst.pdf
