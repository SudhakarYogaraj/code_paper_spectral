.SUFFIXES:

gen_files = tmp/gx.gen tmp/fx.gen tmp/hy.gen tmp/vy.gen tmp/g.gen tmp/f.gen tmp/h.gen tmp/v.gen tmp/lin.gen tmp/rho.gen tmp/fy.gen tmp/drif.gen tmp/diff.gen
sub_files = tmp/gx.sub tmp/fx.sub tmp/hy.sub tmp/vy.sub tmp/g.sub tmp/f.sub tmp/h.sub tmp/v.sub tmp/lin.sub tmp/rho.sub tmp/fy.sub tmp/drif.sub tmp/diff.sub
tmp_files = tmp/gx.tmp tmp/fx.tmp tmp/hy.tmp tmp/vy.tmp tmp/g.tmp tmp/f.tmp tmp/h.tmp tmp/v.tmp tmp/lin.tmp tmp/rho.tmp tmp/fy.tmp tmp/drif.tmp tmp/diff.tmp
target    = ../Problem.cpp

$(target) : problem.init $(tmp_files)
	cat $^ > $@

tmp/%.tmp : aux/%.init tmp/%.sub aux/end
	cat $^ > $@

$(sub_files) : $(gen_files)
	./format_code.sh

$(gen_files) problem.init : build_problem.m
	matlab -nodesktop -nosplash -r build_problem

.PHONY: clean
clean:
	rm -f $(gen_files) $(sub_files) $(tmp_files) $(target)
