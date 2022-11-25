
clean:
	rm -f *.aux *.log *.out *.synctex.gz 
	rm -f exercise??.pdf 

exercises:
	pdflatex exercise00.tex
	pdflatex exercise01.tex
