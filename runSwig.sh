#!/bin/sh

swig -c++ -Isrc -module ChemmineR -o src/r_wrap.cc -r -outdir R swig.i 
