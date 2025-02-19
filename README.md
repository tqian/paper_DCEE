# paper_DCEE

Code for paper "Distal Causal Excursion Effects: Modeling Long-Term Effects of Time-Varying Treatments in Micro-Randomized Trials" by Tianchen Qian

## Files

-   application: code for results in the "Application" section.
-   simulation: code for results in the "Simulation" section.

## How to use the code to replicate results in the paper

-   For results in the "Application" section, follow the steps:
	1. Download data by downloading the entire repository from https://github.com/klasnja/HeartStepsV1 and put it under application/HeartStepsV1-main
	2. Run HeartSteps pre-processing.R
	3. Run HeartSteps_analysis.Rmd

-   For results in the "Simulation" section, follow the steps on a slurm cluster:
	1. Make sure the computing environment is set up (e.g., installing R packages used in the simulation). Modify simulation.slurm as needed (such as the lab account and the R version).
	2. sbatch batch_simulation.sh
	3. Once all simulations are finished (800 jobs, each takes up to 20 minutes on UCI HPC3 cluster), run plot.R to create the plots.


