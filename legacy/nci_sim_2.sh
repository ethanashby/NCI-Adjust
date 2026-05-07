#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eashby@fredhutch.org
ml fhR
Rscript ./nci_sim_2.r