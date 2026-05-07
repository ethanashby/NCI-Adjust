#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eashby@fredhutch.org
#SBATCH --time=30-00:00:00
ml fhR
Rscript ./nci_sim_em.r