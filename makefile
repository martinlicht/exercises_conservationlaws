
default: exercises

clean:
	rm -f *.aux *.log *.out *.synctex.gz 
	rm -f exercise*.pdf 

todolist:
	grep --line-number --color TODO exercise*tex

exercises:
	pdflatex -halt-on-error exercise00_background_and_characteristics.tex
	pdflatex -halt-on-error exercise01.tex
	pdflatex -halt-on-error exercise02.tex
	pdflatex -halt-on-error exercise03.tex
	pdflatex -halt-on-error exercise04.tex
	pdflatex -halt-on-error exercise04a.tex
	pdflatex -halt-on-error exercise05.tex
	pdflatex -halt-on-error exercise06.tex
	pdflatex -halt-on-error exercise06a.tex
	pdflatex -halt-on-error exercise07.tex
	pdflatex -halt-on-error exercise08.tex
	pdflatex -halt-on-error exercise09.tex
	pdflatex -halt-on-error exercise10.tex
	pdflatex -halt-on-error exercise11.tex
	pdflatex -halt-on-error exercise12.tex
	pdflatex -halt-on-error exercise13.tex
