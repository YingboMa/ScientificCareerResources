LL  := latexmk
DEP := $(wildcard *.tex)
MAIN=poster

main: ${DEP}
	${LL} -f -xelatex ${MAIN} -auxdir=output -outdir=output

show:
	zathura output/poster.pdf
