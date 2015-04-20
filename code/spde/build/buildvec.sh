set -v
cd ~/mres/source/program
cat source.m user_input.m > output.m
matlab -nodesktop -nosplash -r output

cd ~/mres/source/program/build
rm built.cpp

# Format result[][] = when not actually a matrix
sed -i 's/^  [^\[]*\[[0-9]\{1,2\}\]\(\[[0-9]\{1,2\}\]\)/    result\1/g' a.gen a_nu.gen drif.gen diff.gen solDrif.gen 

# Format result[][] when the result is actually a matrix
sed -i 's/^  [^\[]*\(\[[0-9]\{1,2\}\]\)\(\[[0-9]\{1,2\}\]\)/    result\1\2/g' dax.gen day.gen solDiff.gen

# When there is only 1 slow process, Matlab produces scalar output. Fix:
sed -i 's/^  [^\[]* =/    result[0] =/g' a.gen a_nu.gen drif.gen diff.gen solDrif.gen
sed -i 's/^  [^\[]* =/    result[0][0] =/g' dax.gen day.gen solDiff.gen

# Matching variable y names
sed -i 's/y\(\([0-9]\|[0-9][0-9]\)\)/y\[\1\]/g' a.gen a_nu.gen dax.gen day.gen drif.gen diff.gen solDrif.gen solDiff.gen  
# In case the slow variable is a vector
sed -i 's/x\(\([0-9]\|[0-9][0-9]\)\)/x\[\1\]/g' a.gen a_nu.gen dax.gen day.gen drif.gen diff.gen solDrif.gen solDiff.gen

cat a_init a.gen end > a.temp
cat a_nu_init a_nu.gen end > a_nu.temp
cat dax_init dax.gen end > dax.temp
cat day_init day.gen end > day.temp
cat drif_init drif.gen end > drif.temp
cat diff_init diff.gen end > diff.temp
cat solDrif_init solDrif.gen end > solDrif.temp
cat solDiff_init solDiff.gen end > solDiff.temp
cat solDrif.temp solDiff.temp a.temp a_nu.temp dax.temp day.temp drif.temp diff.temp > built.cpp

rm *.temp
rm *.gen
