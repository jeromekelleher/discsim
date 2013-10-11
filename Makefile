# simple makefile for development.

SRC=_discsimmodule.c

ext2: ${SRC}
	python setup.py build_ext --inplace

ext3: ${SRC}
	python3 setup.py build_ext --inplace


ctags:
	ctags *.c *.py test/*.py

figs:
	cd docs/asy && make 

docs: ext2 figs 
	cd docs && make clean && make html

	
tags:
	ctags -f TAGS *.c *.py lib/*.[c,h] 
