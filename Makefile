# TEXTARGET is the output filename TEXTSOURCE is input. Simply change these if you want to compile one .tex file 
TEXTARGET=main.pdf
TEXSOURCE=main.tex
# All files on which the TeX-document depends
TEXDEPEND=build/const_a_all.pgf

all: build/$(TEXTARGET)

$(TEXDEPEND): a_plot.py data.txt | build
	python a_plot.py

TEXOPTIONS = -lualatex \
		--output-directory=build \
		--interaction=nonstopmode \
		--halt-on-error \
		--use-make \
		--synctex=1 \

TEXPREFIX = TEXINPUTS=build: \
	BIBINPUTS=build: \
	max_print_line=1048576 \

clean:
	rm -rf build

build:
	mkdir -p build

# Use this to work on tex-document, it will be updated continuesly by latexmk
work: FORCE | build
	$(TEXPREFIX) latexmk -pvc $(TEXOPTIONS) $(TEXSOURCE)

build/$(TEXTARGET): $(TEXSOURCE) $(TEXDEPEND) $(PYTARGET) FORCE | build
	$(TEXPREFIX) latexmk $(TEXOPTIONS) $(TEXSOURCE)

FORCE:

.PHONY: all clean FORCE build
