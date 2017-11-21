# Makefile for the Herring Age Model.

# Targets:
# all - Run everything: 
# 				- compile debug
# 				- compile release
#					- run models id specified directory.
#					- run retrospective analysis in specified directory.
#         - [simulation model test]


all: compile run retro

compile:
	make --directory="src" 

run:
	make --directory="models_2015/sitka" all
	make --directory="models_2015/craig" all

retro:
	make --directory="models_2015/sitka" retro NRET=5
	make --directory="models_2015/craig" retro NRET=5

clean:
	make --directory="src" clean