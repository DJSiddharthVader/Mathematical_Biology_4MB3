To compile the final PDF:

1.
Make sure knitr is your default program to weave Rnw files
(see RStudio > Preferences > Sweave, also select the option ‘pdflatex’ in ‘Typeset LaTeX into PDF using:’)
2.
Open main.Rnw in RStudio
3.
Either press the button “compile PDF” or run the commands:
	library(knitr)
	knit(“main.Rnw”)
4.
If you pressed the button, the PDF should be open in your default viewer. If you ran the commands, open the resulting .tex file with your preferred LaTeX program and compile the PDF.


Notes:

If you want to compile a single question file, include the R chunk at the top of the question file:
<<set-parent, echo=FALSE, cache=FALSE>>=
set_parent('../main.Rnw')
@
