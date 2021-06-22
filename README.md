# Introner Elements

## FIND_CLUSTERS_ALLEUKS.pl

Identifies candidate introner elements in a species by clustering introns based on sequence and length similarity and filtering paralogs. Requires a genome annotation in gff format and corresponding fasta file in the same directory as input. The annotation file and fasta file must have the same prefix (e.g. Name.gff and Name.fa). Outputs candidate introner element familes, separated by family in fasta files ending with ".Pass"

### Dependencies

#### Diamond
Install from https://github.com/bbuchfink/diamond and add to PATH

#### Megablast
Install from https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/ or later and add to PATH


### Usage

```
perl FIND_CLUSTERS_ALLEUKS_7.PL species.gff
```


## cleanTree.py
Removes parenthesis from node names that perturb tree structure.

### Dependencies

**ete3**: http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html

### Usage
```
python cleanTree.py Tree.newick
```
