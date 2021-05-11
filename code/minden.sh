#!/bin/bash

g++ -std=c++11 -o bin/PM code/PM.cpp



#for number in 10 20 50 100 200 500 1000 2000 5000 10000
#do
#bin/PM $number
#done
#exit 0

bin/PM 10000



gnuplot -c 'code/plotter'


