base="$(basename $1 | cut -f1 -d'.')"
if [[ $(which latexmk) ]]; then
    latexmk -c
    R -q -e "library('knitr'); knitr::knit('$1')"
    latexmk -pdf "$base"
    latexmk -pdf  "$base"
else
    R -q -e "library('knitr'); knitr::knit('$1')"
    pdflatex "$base".tex
    biber "$base"
    pdflatex "$base".tex
    pdflatex "$base".tex
fi
