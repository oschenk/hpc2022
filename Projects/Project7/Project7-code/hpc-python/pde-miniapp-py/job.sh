#!/bin/bash
#BSUB -W 00:10
#BSUB -n 4
#BSUB -J pde-miniapp-py

time -p mpirun python main.py 128 100 0.01 v 2>&1 | tee main.log
exit

