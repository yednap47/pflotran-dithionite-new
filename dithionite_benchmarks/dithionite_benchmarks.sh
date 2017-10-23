#!/bin/bash

echo running s2o4-disp
cd s2o4-disp
../../src/pflotran/pflotran -pflotranin s2o4-disp.in >/dev/null 2>/dev/null
python s2o4-disp.py

echo running s2o4-o2
cd ../s2o4-o2
../../src/pflotran/pflotran -pflotranin s2o4-o2.in >/dev/null 2>/dev/null
python s2o4-o2.py

echo running s2o4-fe3
cd ../s2o4-fe3
../../src/pflotran/pflotran -pflotranin s2o4-fe3.in >/dev/null 2>/dev/null
python s2o4-fe3.py

echo running fe2-o2
cd ../fe2-o2
../../src/pflotran/pflotran -pflotranin fe2-o2.in >/dev/null 2>/dev/null
python fe2-o2.py

echo running fe2-cr6
cd ../fe2-cr6
../../src/pflotran/pflotran -pflotranin fe2-cr6.in >/dev/null 2>/dev/null
python fe2-cr6.py

echo 1d_parameterTest
cd ../1d_parameterTest
mpirun -np 8 ../../src/pflotran/pflotran -pflotranin 1d-allReactions-10m-uniformVelocity.in >/dev/null 2>/dev/null
julia plot.jl 1d-allReactions-10m-uniformVelocity

cd ..

if [ -d "results_benchmark" ]; then
  rm -r results_benchmark
fi

mkdir results_benchmark
mv ./*/*.png results_benchmark
echo results are in results_benchmark directory
