SOURCE=M1BIG_bourneuf

all: source clean show

clean:
	- rm *.aux *.log 
	# beamer:
	#- rm *.nav *.snm *.toc *.out *.log *.aux

source:
	pdflatex --enable-write18 $(SOURCE).tex 
	#pdflatex --enable-write18 $(SOURCE).tex && pdflatex $(SOURCE).tex && pdflatex $(SOURCE).tex 


show:
	xdg-open $(SOURCE).pdf


