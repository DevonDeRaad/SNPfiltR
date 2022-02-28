#!/bin/sh
#
#SBATCH --job-name=hippo.qc             # Job Name
#SBATCH --nodes=1             # nodes
#SBATCH --cpus-per-task=1               # CPU allocation per Task
#SBATCH --partition=bi            # Name of the Slurm partition used
#SBATCH --chdir=/home/d669d153/work/dia.din        # Set working d$
#SBATCH --mem-per-cpu=10gb            # memory requested
#SBATCH --time=5000

module load R
R -e "Sys.setenv(RSTUDIO_PANDOC='/panfs/pfs.local/work/bi/bin/pandoc/bin');  rmarkdown::render('qc.Rmd',output_file='qc.html')"
