#!/bin/bash

# script to spawn several instances of matlab,
# in unix, 
unset DISPLAY

for i in 1 2 3 4
do
  #matlab -nodesktop "instance=$i;exdiff3" &
  nohup matlab -r "instance=${i};etdiffnp" >& ${i}_out.log &
  sleep 2
done

