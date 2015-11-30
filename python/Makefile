.SUFFIXES:
.PHONY: clean

# Directories
IN_DIR = inputs
OUT_DIR = ../src/problems
AUX_DIR = tmp

# Files
IN := $(wildcard $(IN_DIR)/*.py)
OUT := $(addprefix $(OUT_DIR)/problem_,$(notdir $(IN:.py=.cpp)))
SUB := $(addprefix $(AUX_DIR)/,$(notdir $(IN:.py=.sub)))
GEN := $(addprefix $(AUX_DIR)/,$(notdir $(IN:.py=.gen)))
BUILD := $(addprefix $(AUX_DIR)/,$(notdir $(IN)))
.SECONDARY: $(OUT) $(SUB)

all : prebuild $(OUT)

prebuild:
	mkdir -p $(AUX_DIR) $(OUT_DIR)

$(OUT_DIR)/problem_%.cpp : problem.init $(AUX_DIR)/%.sub
	cat $^ > $@

$(AUX_DIR)/%.sub : $(AUX_DIR)/%.gen
	sed 's/\<y\(\([0-9]\|[0-9][0-9]\)\)/y\[\1\]/g' $< > $@
	sed -i 's/\<x\(\([0-9]\|[0-9][0-9]\)\)/x\[\1\]/g' $@

$(AUX_DIR)/%.gen : $(AUX_DIR)/%.py
	python $< $@

$(AUX_DIR)/%.py : $(IN_DIR)/%.py
	python read_input.py $< $@

clean :
	rm -f $(SUB) $(SUB) $(GEN) $(BUILD) $(OUT)
