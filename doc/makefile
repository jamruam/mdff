# The compiler
LATEX = pdflatex
BIBTEX = bibtex


EXE= mdff.pdf 

all: $(EXE)

mdff.pdf: 
	$(LATEX) mdff.tex
	$(BIBTEX) mdff.aux
	$(LATEX) mdff.tex
	$(LATEX) mdff.tex        

clean:
	rm -f *.dvi *.pdf
