#!/bin/bash
#$ -l mem_free=10G,h_vmem=10G -cwd -V -j yes -l h_stack=256M

module load conda_R/4.1.x
R CMD BATCH --no-save pop_simulations.R pop_simulations.Rout.${SGE_TASK_ID}

