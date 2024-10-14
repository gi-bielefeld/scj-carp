GENOME = """
## PARAMETERS GENOME LEVEL ###
## Please refer to the manual for more details
#
## Main events
#
## When using the G mode, these events use genome-wise rate.
## The next three events can have gene-family rates when using the Gm mode.
## Rates represent the frequency (and fixation probability) of the different events measure in time units (see your Species Tree)
##
#
DUPLICATION f:{dupl}
TRANSFER f:2
LOSS f:{loss}
#
## The next three events use always genome-wise rates
#
INVERSION f:{inv}
TRANSPOSITION f:{transp}
ORIGINATION f:{orig}
#
## Extensions. Extensions determine how many contiguous genes are affected by a single event
#
DUPLICATION_EXTENSION g:1
TRANSFER_EXTENSION g:1
LOSS_EXTENSION g:1
INVERSION_EXTENSION g:0.05
TRANSPOSITION_EXTENSION g:0.3
#
## Transfer related parameters
#
## REPLACEMENT_TRANSFER gives the probability of the transfer being a replacement transfer. REPLACEMENT_TRANSFER = 0 if all transfers are additive transfers
#
REPLACEMENT_TRANSFER 1
#
# If ASSORTATIVE_TRANSFER is True, the transfer occur preferentially between closely related lineages. The weighted probability is proportional to e ^ - ALPHA * normalized phylogenetic distance.
#  Can be computationally expensive, use at your discretion
#
ASSORTATIVE_TRANSFER False
ALPHA 100
#
## Genome size related parameters
#
INITIAL_GENOME_SIZE	100
#
## MIN_GENOME_SIZE prevents genome to become too small. Any losses affecting genomes under this size will be ignored
#
MIN_GENOME_SIZE	10
#
## Output related parameters. Use SCALE_TREE != 0 only if you used that too in the Species Tree mode. Use the same number
#
EVENTS_PER_BRANCH 1
PROFILES 1
GENE_TREES 1
RECONCILED_TREES 0
VERBOSE 1
SCALE_TREE 0 
#
### GF SPECIFIC PARAMETERS ###
#
GENE_LENGTH f:100
INTERGENE_LENGTH 100
PSEUDOGENIZATION 0.5
#
### GM SPECIFIC PARAMETERS ###
# Include the whole path to the parameters file. The parameters file has a first line (header) and 4 tab separated columns with the name of the gene family (ignore), D, T and L rate.
RATE_FILE False
# If SCALE_RATES = True, it will scale the rates used in the parameters file to the height of the tree, from the leaves to the root node (not the initial node, remember that Zombi also simulates the stem above the root)
SCALE_RATES True
#
## SEED : If 0, the SEED is chosen randomly
SEED 0
"""

TREE = """
### PARAMETERS SPECIES TREE LEVEL ###
## Please refer to the manual for more details
#
SPECIATION	f:1
#
EXTINCTION	f:0.2
#
## STOPPING_RULE = 0 to control depending on TIME, = 1 to control depending on TOTAL_LINEAGES
STOPPING_RULE	1
#
## TOTAL_TIME says when to stop if STOPPING_RULE = 0
#
TOTAL_TIME	1
#
## TOTAL_LINEAGES says when to stop if STOPPING_RULE = 1
#
TOTAL_LINEAGES	20
#
## The simulation fails (and tries again) if the number of lineages in the EST goes below MIN_LINEAGES (when using the stopping rule = 0)
#
MIN_LINEAGES	1
#
## The simulation fails (and tries again) if the number of lineages exceeds MAX_LINEAGES
#
MAX_LINEAGES	10000
#
## VERBOSE = 0 not to print progess, = 1 to do it
#
VERBOSE 0
#
## SCALE_TREE = 0 not scale, SCALE_TREE != 0 scale CROWN to that number
#
SCALE_TREE 0 
#
### Tp SPECIFIC PARAMETERS ###
#
TURNOVER f:0.0002
#
## LINEAGE_PROFILE: TIME1-LINEAGE_GOAL1;TIME2-LINEAGE_GOAL2;...;TIMEk-LINEAGE_GOALk
#
LINEAGE_PROFILE 100-100;300-15000;500-50
#
### Tm SPECIFIC PARAMETERS ###
#
MASS_EXTINCTION 0.8-0.5
#
### Ts SPECIFIC PARAMETERS ###
#
SHIFT_SPECIATION_RATE_FREQUENCY f:1
NUM_SPECIATION_RATE_CATEGORIES 7
BASE_SPECIATION l:1
SHIFT_EXTINCTION_RATE_FREQUENCY f:1
NUM_EXTINCTION_RATE_CATEGORIES 7
BASE_EXTINCTION l:1
#
#
## SEED : If 0, the SEED is chosen randomly
SEED 0
"""

from argparse import ArgumentParser,FileType

def scaled_genome(scale=1):
    dupl=1*scale
    loss=3*scale
    inv=2*scale
    transp=2*scale
    orig=2*scale
    return GENOME.format(dupl=dupl,loss=loss,inv=inv,transp=transp,orig=orig)


parser = ArgumentParser()
parser.add_argument('treeparams',type=FileType('w'))
parser.add_argument('genomeparams',type=FileType('w'))
parser.add_argument('scale',type=float)
args = parser.parse_args()

print(TREE,file=args.treeparams)
print(scaled_genome(scale=args.scale),file=args.treeparams)