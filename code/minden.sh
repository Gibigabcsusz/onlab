#!/bin/bash

g++ -std=c++11 -o bin/PM code/PM.cpp
bin/PM
gnuplot -c 'code/plotter'
