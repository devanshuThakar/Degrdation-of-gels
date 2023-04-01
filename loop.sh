#!/bin/bash

# define time steps to analyze and step size
tstart=$1
tend=$2
tstep=$3

# analysis loop -- loop over LAMMPS time steps
for (( step=$tstart; step<=$tend; step=$step+$tstep ))
do

  echo ${step}
  # write index of largest fragment at timestep step
  ./write_ndx ${step}

  cat afile.txt >> atomfile.txt
  rm afile.txt

done

