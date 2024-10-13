

ZOMBI_HOME = '/prj/gsp/urm/spp_dcj_v2/tools/ZOMBI_01042020/'
ZOMBI_BIN = ZOMBI_HOME+'Zombi.py'
ZOMBI_TREE_PARAMS = ZOMBI_HOME + 'Parameters/SpeciesTreeParameters.tsv'
ZOMBI_GNM_PARAMS = ZOMBI_HOME + 'Parameters/GenomeParameters.tsv'

N_SAMPLES =100
T_RATES = [0.1,0.25,0.5,0.75,1,2.5,5,7.5,10]

rule all_transfer:
    input:
        expand('transferexp/z_transfer{t}/results.txt',t=T_RATES)
    output:
        'transferexp/transfercarp.txt'
    shell:
        'cat {input} > {output}'
        

rule all_default_zombi_genomes:
    input:
        expand('zombi/z_default_sample{s}/G/Genomes/Root_GENOME.tsv',s=range(1,N_SAMPLES+1))

rule all_default_zombi_trees:
    input:
        expand('zombi/z_default_sample{s}/T/ExtantTree.nwk',s=range(1,N_SAMPLES+1))


rule aggregate_treescale_results:
    input:
        #TODO

rule aggregate_transfer_results:
    input:
        expand('carp/measures/z_transfer{{t}}_sample{s}/measure.txt',s=range(1,N_SAMPLES+1))
    output:
        'transferexp/z_transfer{t}/results.txt'
    shell:
        'cat {input} |  sed "s/Carp index:/{wildcards.t}/" > {output}'
    

rule run_carp:
    input:
        pgnm='pangenome/z_{params}/unimog.txt'
    output:
        ci='carp/measures/z_{params}/measure.txt',
        ca='carp/adjacencies/z_{params}/adj.txt'
    log:
        'carp/logs/z_{params}/log.txt'
    shell:
        'python3 scj_carp.py {input} --write-measure {output.ci} --write-carp-adjacencies {output.ca} > {log}'

rule setup_transfer_zombi:
    output:
        tp='zombi/z_transfer{t}_sample{s}/tree_params.tsv',
        gp='zombi/z_transfer{t}_sample{s}/genome_params.tsv'
    shell:
        'cp %s {output.tp};python3 custom_param_files.py {wildcards.t} > {output.gp}'%(ZOMBI_TREE_PARAMS)

rule setup_treescale_zombi:
    output:
        tp='zombi/z_treescale{scl}_sample{s}/tree_params.tsv',
        gp='zombi/z_treescale{scl}_sample{s}/genome_params.tsv'
    shell:
        'python3 treescale.py {output.tp} {output.gp} {wildcards.scl}'

rule setup_default_zombi:
    output:
        tp='zombi/z_default_sample{s}/tree_params.tsv',
        gp='zombi/z_default_sample{s}/genome_params.tsv'
    shell:
        'cp %s {output.tp}; cp %s {output.gp}'%(ZOMBI_TREE_PARAMS,ZOMBI_GNM_PARAMS)
        
        
rule zombi_tree:
    input:
        params='zombi/z_{params}/tree_params.tsv',
    output:
        'zombi/z_{params}/T/ExtantTree.nwk'
    shell:
        'python3 %s T {input.params} zombi/z_{wildcards.params}/'%ZOMBI_BIN


rule zombi_genomes:
    input:
        params='zombi/z_{params}/genome_params.tsv',
        zombitree='zombi/z_{params}/T/ExtantTree.nwk'
    output:
        root='zombi/z_{params}/G/Genomes/Root_GENOME.tsv'
    shell:
        'python3 %s G {input.params} zombi/z_{wildcards.params}/'%ZOMBI_BIN


rule zombi_to_pangenome:
    input:
        rgnm='zombi/z_{params}/G/Genomes/Root_GENOME.tsv',
        tree='zombi/z_{params}/T/ExtantTree.nwk'
    output:
        'pangenome/z_{params}/unimog.txt'
    shell:
        'python3 zombi2unimog.py {input.tree} zombi/z_{wildcards.params}/G/Genomes/ > {output}'
