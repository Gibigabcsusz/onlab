#!/bin/bash

g++ -std=c++11 -o bin/PM code/PM.cpp
bin/PM
gnuplot -e "set terminal eps; set output 'output/plot.eps'; plot 'output/semmi.dat' with lines"
