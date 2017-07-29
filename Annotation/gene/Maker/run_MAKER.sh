#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=40G
#SBATCH --time=9-00:00:00
#SBATCH --output=fc-maker-%A-%a.out
#SBATCH --job-name="fc-maker"
#SBATCH --array=1-20
#SBATCH -p intel

module unload perl
module load maker/2.31.8

maker
