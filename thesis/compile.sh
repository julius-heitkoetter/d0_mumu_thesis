# make the output directory and copy everything into it
mkdir -p output
cp * output
cp -r fontsets output
cd output

# Compile the pdf document the first time
pdflatex  thesis.tex

# Compile the bibliography
bibtex thesis

# Compile the pdf document the second time
pdflatex thesis.tex

# Compile the pdf document the third time...
pdflatex thesis.tex

cd ..