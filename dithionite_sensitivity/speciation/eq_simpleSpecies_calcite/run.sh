#!/bin/bash

../../../src/pflotran/pflotran -pflotranin calcite.in
../../../src/pflotran/pflotran -pflotranin calcite_initial.in
../../../src/pflotran/pflotran -pflotranin calcite_injectant.in

