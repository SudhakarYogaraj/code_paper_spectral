.PHONY : all clean

all := summary.pdf

summary.pdf : summary.yaml
	python convert.py summary.yaml summary.tex
	latexmk summary.tex

clean : 
	latexmk -c
	rm summary.pdf
	rm summary.tex
