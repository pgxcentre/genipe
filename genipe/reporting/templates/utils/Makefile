PROJECT=report
TEX=pdflatex
BIBTEX=bibtex
BUILDTEX=$(TEX) $(PROJECT).tex

all:
	$(BUILDTEX)
	$(BIBTEX) $(PROJECT)
	$(BUILDTEX)
	$(BUILDTEX)

clean-all:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc $(PROJECT).pdf

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc
