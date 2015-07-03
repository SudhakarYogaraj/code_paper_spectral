.SUFFIXES:
.PHONY: clean

target = ../src/Problem.cpp
sub = tmp/built.cpp
gen = tmp/output.gen

$(target) : problem.init $(sub)
	cat $^ > $@

$(sub) : $(gen)
	sed 's/\<y\(\([0-9]\|[0-9][0-9]\)\)/y\[\1\]/g' $(gen) > $(sub)
	sed -i 's/\<x\(\([0-9]\|[0-9][0-9]\)\)/x\[\1\]/g' $(sub)

$(gen) : build_problem.py
	python $<

clean:
	rm -f $(gen) $(sub) $(target)
