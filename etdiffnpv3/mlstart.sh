#!/bin/bash

unset DISPLAY

nohup matlab -nodesktop -nojvm -r etmulti >& et1.log &

