.SUFFIXES:
.PHONY: clean
.PRECIOUS: *.tex

# Directories
IN_DIR = inputs
OUT_DIR = ../src/problems
AUX_DIR = tmp
TEX_DIR = tex

# Files
IN := $(wildcard $(IN_DIR)/*.py)
OUT := $(addprefix $(OUT_DIR)/,$(notdir $(IN:.py=.cpp)))
SUB := $(addprefix $(AUX_DIR)/,$(notdir $(IN:.py=.sub)))
GEN := $(addprefix $(AUX_DIR)/,$(notdir $(IN:.py=.gen)))
TEX := $(addprefix $(TEX_DIR)/,$(notdir $(IN:.py=.tex)))
BUILD := $(addprefix $(AUX_DIR)/,$(notdir $(IN)))
.SECONDARY: $(OUT) $(SUB)

all : prebuild $(OUT)

prebuild:
	mkdir -p $(AUX_DIR) $(OUT_DIR) $(TEX_DIR)

$(OUT_DIR)/%.cpp : problem.init $(AUX_DIR)/%.sub
	cat $^ > $@

$(AUX_DIR)/%.sub : $(AUX_DIR)/%.gen
	sed 's/\<y\(\([0-9]\|[0-9][0-9]\)\)/y\[\1\]/g' $< > $@
	sed -i 's/\<x\(\([0-9]\|[0-9][0-9]\)\)/x\[\1\]/g' $@

$(AUX_DIR)/%.gen: $(AUX_DIR)/%.py
	python $< $@ $(TEX_DIR)/$(notdir $(@:.gen=.tex))

$(AUX_DIR)/%.py : $(IN_DIR)/%.py
	python read_input.py $< $@

clean :
	rm -f $(SUB) $(SUB) $(GEN) $(BUILD) $(OUT) $(TEX)
