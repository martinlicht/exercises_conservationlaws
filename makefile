
default: exercises

clean:
	rm -f *.aux *.log *.out *.synctex.gz 
	rm -f exercise??.pdf 

exercises:
	pdflatex exercise00_background_and_characteristics.tex
	pdflatex exercise01.tex
	pdflatex exercise02.tex
	pdflatex exercise03.tex
	pdflatex exercise04.tex
