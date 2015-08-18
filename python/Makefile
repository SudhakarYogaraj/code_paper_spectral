.SUFFIXES:
.PHONY: clean

target = ../cpp/src/Problem.cpp
sub = tmp/built.cpp
gen = tmp/output.gen
build = build_problem.py
input = user_input.py build_problem_init.py

all :
	mkdir -p tmp
	make target

target : $(target)

$(target) : problem.init $(sub)
	cat $^ > $@

$(sub) : $(gen)
	sed 's/\<y\(\([0-9]\|[0-9][0-9]\)\)/y\[\1\]/g' $(gen) > $(sub)
	sed -i 's/\<x\(\([0-9]\|[0-9][0-9]\)\)/x\[\1\]/g' $(sub)

$(gen) : $(build)
	python $<

$(build) : $(input)
	python read_input.py

clean:
	rm -f $(gen) $(sub) $(target)
