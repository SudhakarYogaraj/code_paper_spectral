#!/bin/bash
cd tmp

for file in *.gen; do
    cp $file `echo $file | sed 's/gen/sub/g'`
done

dou="v.sub lin.sub rho.sub"
vec="vy.sub g.sub f.sub h.sub drif.sub diff.sub"
mat="gx.sub fx.sub hy.sub fy.sub"

# Formatting for double, vector, and matrix outputs
sed -i 's/^  [^=]*/    result /g' $dou
sed -i 's/^  [^\[]*\[[0-9]\{1,2\}\]\(\[[0-9]\{1,2\}\]\)/    result\1/g' $vec
sed -i 's/^  [^\[]*\(\[[0-9]\{1,2\}\]\)\(\[[0-9]\{1,2\}\]\)/    result\1\2/g' $mat

# Fix of the case where there is only 1 slow process
sed -i 's/^  [^\[]* =/    result[0] =/g' $vec
sed -i 's/^  [^\[]* =/    result[0][0] =/g' $mat

# Matching c++ program variables
sed -i 's/y\(\([0-9]\|[0-9][0-9]\)\)/y\[\1\]/g' *.sub
sed -i 's/x\(\([0-9]\|[0-9][0-9]\)\)/x\[\1\]/g' *.sub
