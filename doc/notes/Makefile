all:
	pandoc -s --bibliography cite.bib --toc --template=template.latex\
		--filter pandoc-eqnos \
		--filter pandoc-citeproc \
		--number-sections text.md -o text.pdf #--latex-engine=xelatex 
# If Chinese language is needed, uncomment the above and 
# add this line in .latex file 
