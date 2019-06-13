# Makefile for compiling the model.tpl code
EXEC  = ham
TPL  := $(EXEC).tpl
DAT  := sitka_2.ctl sitka_2.dat
NSIMS = 10






$(EXEC): $(TPL)
	@admb -g $(EXEC)

$(EXEC).par: $(EXEC) 
	./$(EXEC) $(ARG)

.PHONY: run runs
run: $(EXEC).par

runmc:
	./$(EXEC) $(ARG)




simdirs := $(shell echo 'cat(formatC(1:$(NSIMS),digits=3,flag="0"))' | R --vanilla --slave)
setsims := $(foreach dir,$(simdirs),$(dir)/setsims)
runsims := $(foreach dir,$(simdirs),$(dir)/runsims)

$(setsims): $(EXEC)
	mkdir $(@D);
	cp $(EXEC) $(EXEC).dat $(DAT) Makefile $(@D)
	cd $(@D); touch setsims

set: $(setsims)


$(runsims):
	cd $(@D); make runmc ARG='-sim -$(@D)'

sim: set $(runsims)
	touch sim


clean:
	rm -r 0* \
	$(EXEC) $(EXEC).[bcehlprs]* *[0-9].cpp *[0-9].htp *.obj \
	admb* admodel* fmin.log variance *.prj

# Develop an R-script to read and summarise the data from simulations.

READSIMS = 'source(file.path("../R/readSims.R"))'
RENDEREP = 'library(rmarkdown); render("../R/ham.Rmd",output_file="ham2.html")'

sims.Rdata: sim
	@echo $(READSIMS) | R --vanilla --slave
	@echo $(RENDEREP) | R --vanilla --slave

collect: sims.Rdata

## end of makefile