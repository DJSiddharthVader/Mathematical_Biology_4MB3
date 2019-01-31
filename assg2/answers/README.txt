To compile the final PDF:

1.
Make sure knitr is your default program to weave Rnw files 
(see RStudio > Preferences > Sweave, uncheck ‘Always enable Rnw concordance’)

2. 
Open math4MB3_A2_PlagueDoctors.Rnw in RStudio

3. 
Either press the button “compile PDF” or run the commands:
	library(knitr)
	knit(“math4MB3_A2_PlagueDoctors.Rnw”)

4. 
If you pressed the button, the PDF should be open in your default viewer. If you ran the commands, open the resulting .tex file with your preferred LaTeX program and compile the PDF.

	
Notes:

For some reason, the question files with plots don’t compile well if they are in a subfolder, so they must be in the same directory as the main document (errors occur when dev=‘tikz’ is enabled in an R chunk within a subfolder/child-file.Rnw).

If you want to compile a single question file, include the R chunk at the top of the question file:
<<set-parent, echo=FALSE, cache=FALSE>>=
set_parent('math4MB3_A2_PlagueDoctors.Rnw')
@