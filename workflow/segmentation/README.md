# Segmentation Workflow
=====================================

## Introduction
---------------

This workflow evaluates the SCJ-CARP measure under different block definitions on several pangenomes. The original setup is described in the Section "Evaluating the SCP-CARP Measure of Multiple Pangenomes" in the manuscript.

## Configuration
---------------

To use this workflow, you need to adapt the `config.yaml` file to your needs. The config file should contain the following parameters:
* `run_list`: list of runs to process
* `runs`: dictionary of run configurations
* `k_sibeliaz`: k-mer size for Sibeliaz (default: 25)
* `a_sibeliaz`: list of alpha values for Sibeliaz (default: [2*len(genomeFiles), 2*10*len(genomeFiles)])

## Dependencies
---------------

This workflow depends on [Minigraph-Cactus](https://doi.org/10.1038/s41587-023-01793-w), [Sibeliaz](https://doi.org/10.1038/s41467-020-19777-8) and [maf2synteny](https://doi.org/10.1101/gr.236273.118). Make sure they are installed and adapt the wrapper scripts `adapt_*.sh` to your installation. Then remove the `adapt_` prefix.

## Data Preparation
---------------

Download your data and set up your `config.yaml`. An example is given in `config.yaml` and the original configuration in `original.yaml`. You can use the script `autonames.py` to help with data preparation. Make sure all fasta-files contain only `ACGTN`-characters, not IUPAC as this causes errors in `cactus` otherwise.

## Running the Workflow
-------------------

The workflow can then be run with `snakemake --cores <num-your-available-cores>`.

## Post-processing
-----------------

`collect_results.py`, `plot_results.py`, and `plot_results_2d.py` may help you evaluate the results.

## Workflow Overview
-------------------

The target rule `all_carp` generates the CARP measure files for all specified runs and tools. The results can be found in the `carp` directory within the specified output directory.
