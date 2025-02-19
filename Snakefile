

ZOMBI_HOME = '/prj/gsp/urm/spp_dcj_v2/tools/ZOMBI_01042020/'
ZOMBI_BIN = ZOMBI_HOME+'Zombi.py'
ZOMBI_TREE_PARAMS = ZOMBI_HOME + 'Parameters/SpeciesTreeParameters.tsv'
ZOMBI_GNM_PARAMS = ZOMBI_HOME + 'Parameters/GenomeParameters.tsv'

N_SAMPLES = 100
T_RATES = [0.125,0.25,0.5,1,2,4,8,16,32]
TREE_SCALES = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]


rule all_treescale_results:
    input:
        expand('treescale/z_treescale{scl}_sample{s}/agg.txt',scl=TREE_SCALES,s=range(1,N_SAMPLES+1))
    output:
        'treescale/treescale.txt'
    shell:
        'cat {input} > {output}'

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


rule cargo_build:
    output:
        './target/release/scj-carp-rust'
    shell:
        'cargo build -r'


rule aggregate_treescale_results_hor:
    input:
        pr='pr/z_treescale{scl}_sample{s}/pr.txt',
        ci='carp/measures/z_treescale{scl}_sample{s}/measure.txt'
    output:
        'treescale/z_treescale{scl}_sample{s}/agg.txt'
    shell:
        'echo {wildcards.scl},$(cat {input.pr}),$(cat {input.ci} | sed "s/Carp index: //") > {output}'

rule aggregate_transfer_results:
    input:
        expand('transferexp/z_transfer{{t}}_sample{s}/summary.txt',s=range(1,N_SAMPLES+1))
    output:
        'transferexp/z_transfer{t}/results.txt'
    shell:
        'cat {input} > {output}'

rule transfer_common:
    input:
        ms='carp/measures/z_transfer{t}_sample{s}/measure.txt',
        pr='pr/z_transfer{t}_sample{s}/pr.txt'
    output:
        'transferexp/z_transfer{t}_sample{s}/summary.txt'
    shell:
        'echo $(cat {input.ms} |  sed "s/Carp index: /{wildcards.t},/"),$(cat {input.pr}) > {output}'
    

rule prec_recall:
    input:
        ca='carp/adjacencies/z_{params}/adj.txt',
        gt='pangenome/z_{params}/adj_root.txt'
    output:
        'pr/z_{params}/pr.txt'
    shell:
        'python3 precrecall.py --ground {input.gt} --result {input.ca} > {output}'

rule run_carp:
    input:
        pgnm='pangenome/z_{params}/unimog.txt',
        binary = 'target/release/scj-carp-rust'
    output:
        ci='carp/measures/z_{params}/measure.txt',
        ca='carp/adjacencies/z_{params}/adj.txt'
    log:
        'carp/logs/z_{params}/log.txt'
    shell:
        '{input.binary} --unimog {input.pgnm} --write-measure {output.ci} --write-ancestor {output.ca} > {log}'

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
    log:
        'zombi/z_{params}/T/log.txt'
    shell:
        'python3 %s T {input.params} zombi/z_{wildcards.params}/ > {log} 2>&1 '%ZOMBI_BIN


rule zombi_genomes:
    input:
        params='zombi/z_{params}/genome_params.tsv',
        zombitree='zombi/z_{params}/T/ExtantTree.nwk'
    output:
        root='zombi/z_{params}/G/Genomes/Root_GENOME.tsv'
    log:
        'zombi/z_{params}/G/log.txt'
    shell:
        'python3 %s G {input.params} zombi/z_{wildcards.params}/ > {log} 2>&1 '%ZOMBI_BIN


rule zombi_to_pangenome:
    input:
        rgnm='zombi/z_{params}/G/Genomes/Root_GENOME.tsv',
        tree='zombi/z_{params}/T/ExtantTree.nwk'
    output:
        'pangenome/z_{params}/unimog.txt'
    shell:
        'python3 zombi2unimog.py {input.tree} zombi/z_{wildcards.params}/G/Genomes/ > {output}'


rule root_adjacencies:
    input:
        'pangenome/z_{params}/root.txt'
    output:
        'pangenome/z_{params}/adj_root.txt'
    shell:
        'python3 genome2adj.py {input} > {output}'


rule zombi_to_root:
    input:
        rgnm='zombi/z_{params}/G/Genomes/Root_GENOME.tsv',
        tree='onlyroot.nwk'
    output:
        'pangenome/z_{params}/root.txt'
    shell:
        'python3 zombi2unimog.py {input.tree} zombi/z_{wildcards.params}/G/Genomes/ --onlyRoot > {output}'
