#!/bin/bash

no_perms=1000
echo $no_perms " rotated networks to be processed"

qsub -l mem_free=0.5G -binding linear:1 -N rot-classifictaion -t 1-${no_perms} -j y ./q_from_array.sh

