GENOME_TEMPLATE_TRANSFER="""### PARAMETERS GENOME LEVEL ###
## Please refer to the manual for more details
#
## Main events
#
## When using the G mode, these events use genome-wise rate.
## The next three events can have gene-family rates when using the Gm mode.
## Rates represent the frequency (and fixation probability) of the different events measure in time units (see your Species Tree)
##
#
DUPLICATION f:1
TRANSFER f:{transfer}
LOSS f:3
#
## The next three events use always genome-wise rates
#
INVERSION f:2
TRANSPOSITION f:2
ORIGINATION f:2
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
INITIAL_GENOME_SIZE	1000
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

from argparse import ArgumentParser
parser=ArgumentParser()
parser.add_argument('transfer',type=float)

args = parser.parse_args()

print(GENOME_TEMPLATE_TRANSFER.format(transfer=args.transfer))
