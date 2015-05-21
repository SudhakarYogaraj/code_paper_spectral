#!/bin/bash

# Directory of the code
cd ~/Dropbox/phd/papers/spectral/code/finite

# Launch matlab program
matlab -nodesktop -nosplash -r build_problem_1d

# Directory of tmp files
cd tmp/

# Deletion of useless files
rm -f *.temp
rm -f built.cpp

# Substitution of variable names
sed -i 's/\<x\>/x[0]/g' fx.gen hy.gen vy.gen g.gen f.gen h.gen v.gen lin.gen rho.gen
sed -i 's/\<y\>/y[0]/g' fx.gen hy.gen vy.gen g.gen f.gen h.gen v.gen lin.gen rho.gen

sed -i 's/^  [^=]*/    result[0][0] /g' fx.gen hy.gen
sed -i 's/^  [^=]*/    result[0] /g' vy.gen g.gen f.gen h.gen
sed -i 's/^  [^=]*/    result /g' v.gen lin.gen rho.gen

cat fx.init fx.gen end > fx.temp
cat hy.init hy.gen end > hy.temp
cat vy.init vy.gen end > vy.temp
cat g.init g.gen end > g.temp
cat f.init f.gen end > f.temp
cat h.init h.gen end > h.temp
cat v.init v.gen end > v.temp
cat lin.init lin.gen end > lin.temp
cat rho.init rho.gen end > rho.temp

cat fx.temp hy.temp vy.temp g.temp f.temp h.temp v.temp lin.temp rho.temp > built.cpp
cat ../aux_problem.cpp built.cpp > ../Problem.cpp
