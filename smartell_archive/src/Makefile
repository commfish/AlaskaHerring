LIB:=
INC:=
OUT:=ham

ifdef ComSpec
	RM=del /F /Q
	COPY=copy 
else
	RM=rm -rf
	COPY=cp
endif


all: build/debug/ham build/release/ham

build/debug/ham: build/debug $(OUT).tpl $(LIB) $(INC)
	@admb -g $(OUT).tpl
	@$(COPY) $(OUT) build/debug/
	@$(RM) $(OUT)

build/debug:
	@echo creating debug directory
	@mkdir -p build/debug


build/release/ham: build/release $(OUT).tpl $(LIB) $(INC)
	@admb -g $(OUT).tpl
	@$(COPY) $(OUT) build/release/
	@$(RM) $(OUT)

build/release:
	@echo creating release directory
	@mkdir -p build/release


clean:
	@$(RM) build $(OUT).cpp $(OUT).htp $(OUT).o*