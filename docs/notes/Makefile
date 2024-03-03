# The document, distributions.tex, relies on the jlstyle Latex package in order
# to show Julia code inline. You can install that package into your latex,
# or it is also included here in the subdirectory called texmf. Setting the
# TEXINPUTS variable is one way to tell pdflatex where to find TEXINPUTS.
distributions.pdf: distributions.tex
	TEXINPUTS=./texmf//:${TEXINPUTS} pdflatex distributions
	TEXINPUTS=./texmf//:${TEXINPUTS} pdflatex distributions
