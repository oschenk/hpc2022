filename=template

pdf:
	pdflatex ${filename}
	pdflatex ${filename}
	make clean

read:
	evince ${filename}.pdf &


clean:
	rm -f ${filename}.out ${filename}.log ${filename}.bbl ${filename}.blg ${filename}.au ${filename}.log  ${filename}.ps ${filename}.aux ${filename}.out ${filename}.dvi ${filename}.bbl ${filename}.blg ${filename}.toc  ${filename}.nav ${filename}.vrb ${filename}.snm
