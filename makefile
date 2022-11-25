
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
	pdflatex exercise04a.tex
	pdflatex exercise05.tex
	pdflatex exercise06.tex
	pdflatex exercise06a.tex
	pdflatex exercise07.tex
	pdflatex exercise08.tex
	pdflatex exercise09.tex
	pdflatex exercise10.tex
