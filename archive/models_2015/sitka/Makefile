EXEC=ham
DAT=sitka_2.dat
CTL=sitka_2.ctl
PIN=sim0.pin

TPLS=../../src/$(EXEC).tpl
DATA=$(EXEC).dat
SRCS=../../src/build
ARGS=


ifdef DEBUG
	DIST=$(SRCS)/debug/$(EXEC)
else
	DIST=$(SRCS)/release/$(EXEC)
endif


all: $(TPLS) $(EXEC) $(EXEC).rep

$(EXEC): $(DIST)
	ln -s $(DIST) $@


$(EXEC).rep: $(EXEC) $(DATA) $(DAT) $(CTL)
	./$(EXEC) -ind $(DATA) $(ARGS)

$(TPLS):
	make --directory "../../"

run: $(DATA) 
	./$(EXEC) -ind $(DATA) $(ARGS)


mcmc: $(EXEC) $(DATA) $(CTL) $(DAT)
	./$(EXEC) -ind $(DATA) $(ARGS) -mcmc -mcsave 100
	./$(EXEC) -ind $(DATA) $(ARGS) -mceval

# Retrospective Rule
NRET=1
retdirs := $(shell echo 'cat(formatC(1:$(NRET),digits=3,flag="0"))' | R --vanilla --slave)
makdirs := $(foreach dir, $(retdirs), $(dir)/makdirs)
rundirs := $(foreach dir, $(retdirs), $(dir)/rundirs)

$(makdirs): 
	mkdir -p rets/$(@D)
	cp $(EXEC) $(DATA) $(DAT) $(CTL) Makefile rets/$(@D);
	cd $(@D); touch makdirs

$(rundirs): 
	@cd rets/$(@D); make run ARGS='-retro $(@D)'
	@cd rets/$(@D); touch rundirs

retro: $(EXEC) $(makdirs) $(rundirs)
	touch RETRODONE

# Simulation Rule
NSIMS=2
simdirs := $(shell echo 'cat(formatC(1:$(NSIMS),digits=3,flag="0"))' | R --vanilla --slave)
setsims := $(foreach dir,$(simdirs),$(dir)/setsims)
runsims := $(foreach dir,$(simdirs),$(dir)/runsims)

$(setsims):
	mkdir -p sims/$(@D);
	cp $(EXEC) $(DATA) $(DAT) $(CTL) $(PIN) Makefile sims/$(@D);
	cd $(@D); touch setsims

$(runsims):
	cd sims/$(@D); make run ARGS="-sim -$(@D) -ainp $(PIN)";

sim: $(EXEC) $(setsims) $(runsims)
	touch SIMSDONE

dust:
	@rm -f admodel.cov
	@rm -f admodel.hes
	@rm -f admodel.dep
	@rm -f fmin.log variance

clean: dust
	@rm -f $(EXEC) 
	@rm -f $(EXEC).rep $(EXEC).par ham.[bcehmrpls]*
	@rm -r rets ||:
	@rm -r sims ||:
	@rm -f makedirs
	@rm -f retro
	@rm -f RETRODONE
	@rm -f SIMSDONE
	@rm -f ders.dat makdirs derivatives 
	@rm -f ssb.ps
