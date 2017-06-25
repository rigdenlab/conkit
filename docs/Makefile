# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS       =
SPHINXBUILD      = sphinx-build
PAPER            =
BUILDDIR         = _build
GH_PAGES_SOURCES = conkit docs

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don't have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) .

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest gettext gh-pages

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  epub       to make an epub"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"
	@echo "  figures    to create all figures used in the documentation"
	@echo "  apidoc     to create the apidoc pages"

clean:
	rm -rf $(BUILDDIR)/*

html:
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

singlehtml:
	$(SPHINXBUILD) -b singlehtml $(ALLSPHINXOPTS) $(BUILDDIR)/singlehtml
	@echo
	@echo "Build finished. The HTML page is in $(BUILDDIR)/singlehtml."

epub:
	$(SPHINXBUILD) -b epub $(ALLSPHINXOPTS) $(BUILDDIR)/epub
	@echo
	@echo "Build finished. The epub file is in $(BUILDDIR)/epub."

doctest:
	$(SPHINXBUILD) -b doctest $(ALLSPHINXOPTS) $(BUILDDIR)/doctest
	@echo "Testing of doctests in the sources finished, look at the " \
	      "results in $(BUILDDIR)/doctest/output.txt."

figures:
	if [ -d _tmp ]; then rm -rf _tmp; else mkdir _tmp; fi
	git clone https://github.com/rigdenlab/conkit-examples.git _tmp/conkit-examples 
	conkit-plot chord --confidence -o examples/images/toxd_chord_simple.png _tmp/conkit-examples/toxd/toxd.fasta fasta _tmp/conkit-examples/toxd/toxd.mat ccmpred
	conkit-plot cmap -o examples/images/toxd_cmap_simple.png _tmp/conkit-examples/toxd/toxd.fasta fasta _tmp/conkit-examples/toxd/toxd.mat ccmpred
	conkit-plot cmap -o examples/images/toxd_cmap_reference.png -p _tmp/conkit-examples/toxd/toxd.pdb -pf pdb _tmp/conkit-examples/toxd/toxd.fasta fasta _tmp/conkit-examples/toxd/toxd.mat ccmpred
	conkit-plot cmap -o examples/images/toxd_cmap_advanced.png -e _tmp/conkit-examples/toxd/toxd.psicov -ef psicov -p _tmp/conkit-examples/toxd/toxd.pdb -pf pdb _tmp/conkit-examples/toxd/toxd.fasta fasta _tmp/conkit-examples/toxd/toxd.mat ccmpred
	conkit-plot cmap -o examples/images/toxd_cmap_confidence.png --confidence -e _tmp/conkit-examples/toxd/toxd.psicov -ef psicov -p _tmp/conkit-examples/toxd/toxd.pdb -pf pdb _tmp/conkit-examples/toxd/toxd.fasta fasta _tmp/conkit-examples/toxd/toxd.mat ccmpred
	conkit-plot scov -o examples/images/toxd_scov_plot.png _tmp/conkit-examples/toxd/toxd.a3m a3m
	conkit-plot peval -o examples/images/toxd_peval_plot.png -j 0.1 -min 0 -max 2 _tmp/conkit-examples/toxd/toxd.pdb pdb _tmp/conkit-examples/toxd/toxd.fasta fasta _tmp/conkit-examples/toxd/toxd.mat ccmpred
	conkit-plot cdens -o examples/images/4p9g_cdens_plot.png _tmp/conkit-examples/4p9g/4p9g.fasta fasta _tmp/conkit-examples/4p9g/4p9g.mat ccmpred
	rm -rf _tmp

apidoc:
	sphinx-apidoc -e -f -l -M -o api/generated ../conkit
	
