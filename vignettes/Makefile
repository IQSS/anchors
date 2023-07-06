# Makefile for LaTeX files

LATEX	= pdflatex
BIBTEX	= bibtex
DVIPS	= dvips
PSTOPDF = ps2pdf

SRC	= anchors.tex

all 	: anchors.pdf

%.pdf:	%.tex
	$(LATEX) $*
	$(BIBTEX) $*
	$(LATEX) $*
	$(LATEX) $*
#	$(DVIPS) $*
#	$(PSTOPDF) $*.ps
#	rm -f *.ps anchors-* $*.dvi $*.log $*.bbl *.tex $*.aux $*.blg $*.log
	rm -f $*.dvi $*.log $*.aux $*.blg  $*.out $*.bbl 

