# Simulated Experiments Workflow


## Introduction


This workflow performs simulated experiments to evaluate the performance of the CARP measure under different scenarios. The original setup is described in the Section "Simulated Experiments" in the manuscript.

## Configuration


To use this workflow, you need to adapt the `config.yaml` file to your needs. The config file should contain the following parameters:
* `zombi_home`: path to the Zombi installation
* `n_samples`: number of samples to generate per step
* `transfer_rates`: list of transfer rates to simulate
* `tree_scales`: list of rearrangement rates to simulate

## Dependencies


This workflow depends on [Zombi](https://doi.org/10.1093/bioinformatics/btz710).

## Running the Workflow


The workflow can be run with `snakemake <target-rule> --cores <num-your-available-cores>`.

## Post-processing


The resulting reconstruction precision and recall can be plotted using the scripts `rearrangement_vs_precrec.py` and `transfer_vs_precrec.py` in the `plotscripts` directory at the root of the repository. The correlation between the measure and operations can be plotted using the scripts `rearrangement_vs_measure.py` and `transfer_vs_measure.py`.

## Workflow Overview


The target rules `all_treescale` and `all_transfer` generate the results for all specified tree scales and transfer rates. The results can be found in the `treescale` and `transferexp` directories, respectively.
