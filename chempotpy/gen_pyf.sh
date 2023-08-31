#!/bin/bash

# Script to generate .pyf file for new surfaces
# Dayou Zhang, Aug 30, 2023
 
src=`basename $1`
module=${src%.*}
cd `dirname $1`
f2py only: pes : -m ${module} -h ${module}.pyf $src
