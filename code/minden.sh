#!/bin/bash

g++ -std=c++11 -o bin/PM code/PM.cpp

touch output/meres_ng
> output/meres_ng

for number in 10 20 50 100 200 500 1000 2000 5000
do
bin/PM $number 300 >> output/meres_ng
done
exit 0




#gnuplot -c 'code/plotter'

gnuplot set terminal eps; set output 'output/meres_ng.eps'; plot 'output/meres_ng' with lines gnuplot set terminal eps;
