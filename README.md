# Introner Elements

## Tools
#### FIND_CLUSTERS_ALLEUKS_7.PL

Identifies candidate introner elements in a species by clustering introns based on sequence and length similarity and filtering paralogs. Requires a genome annotation in gff format and corresponding fasta file as input. The annotation file and fasta file must have the same prefix (e.g. Name.gff and Name.fa).

```
perl FIND_CLUSTERS_ALLEUKS_7.PL species.gff
```

## Scripts

#### cleanTree.py

Removes parenthesis from node names that perturb tree structure.

```
python cleanTree.py Tree.newick
```
