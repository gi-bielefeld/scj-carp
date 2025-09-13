# Segmentation Workflow

This worklfow evaluates the SCJ-CARP measure under different block definitions on several pangenomes. The original setup is described in the Section "Evaluating the SCP-CARP Measure of Multiple Pangenomes" in the manuscript.

This workflow depends on cactus and sibeliaz. Make sure they are installed and adapt the wrapper scripts
`adapt_*.sh` to your installation. Then remove the `adapt_` prefix.

Then download your data and set up your `config.yaml`. An example is given in `config.yaml` and the original configuration in `original.yaml`

You can use the script `autonames.py`. Make sure all fasta-files contain only `ACGTN`-characters, not IUPAC as this causes errors in cactus otherwise.

The workflow can then be run with `snakemake --cores <num-your-available-cores>`.

`collect_results.py`,`plot_results.py` and `plot_results_2d.py` may help you evaluate the results.