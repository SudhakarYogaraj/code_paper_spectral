#!/bin/bash
cd tmp

for file in *.gen; do
    cp $file `echo $file | sed 's/gen/sub/g'`
done

# Substitution of variable names
sed -i 's/\<x\>/x[0]/g' fx.sub hy.sub vy.sub g.sub f.sub h.sub v.sub lin.sub rho.sub
sed -i 's/\<y\>/y[0]/g' fx.sub hy.sub vy.sub g.sub f.sub h.sub v.sub lin.sub rho.sub

# Formating of output data type
sed -i 's/^  [^=]*/    result[0][0] /g' fx.sub hy.sub
sed -i 's/^  [^=]*/    result[0] /g' vy.sub g.sub f.sub h.sub
sed -i 's/^  [^=]*/    result /g' v.sub lin.sub rho.sub
