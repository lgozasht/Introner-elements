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

## GetCoordinates.pl
Appends coordinates for IEs identified by "FIND_CLUSTERS_ALLEUKS.pl". Run in the same directory as the inpput and output files of "FIND_CLUSTERS_ALLEUKS.pl". Produces output files for each family ending in ".withcoords".

### Usage

```
perl GetCoordinates.pl *.Pass
```

## blastBack.py
Constructs a consensus sequence for each IE families using a positional weight matrix and blasts each consensus back at the reference to identify candidate IEs that exist outside of genes or are not annotated as introns. Also converts IE family file names into numbered fasta files. Requires a tab delimitted file containing the directory path and reference fasta file (in that order) for each considered genome separated by line. Each directory should contain all IE family ".withcoords" files for a given genome and the respective reference fasta file.

### Dependencies

#### BLAST+
Install from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ and add to PATH

#### MAFFT
Install from https://mafft.cbrc.jp/alignment/software/source.html and add to PATH

### Usage

```
python blastBack.py directory_list.tsv
```

## outsideGenes.py
Uses output from "blastBack.py" to identify candidate IE insertions outside of genes and unannotated insertions inside of genes. Deposits these IEs in separate fasta files. Requires a tab delimitted file containing the directory path, the reference fasta file and the annotation file (in that order) for each considered genome separated by line. Each directory should contain all IE family files for a given genome, the output from "blastBack.py" and the respective reference fasta file and annotation.

### Usage

```
python outsideGenes.py directory_list.tsv
```

## removeDuplicates.py
Removes duplicate hits across IE families after the blast back step. Requires a text file with the path to each directory containing the output of "outsideGenes.py"  separated by line.

### Usage

```
python removeDuplicates.py directory_list.tsv
```


## cleanTree.py
Removes parenthesis from node names that perturb tree structure.

### Dependencies

**ete3**: http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html

### Usage
```
python cleanTree.py Tree.newick
```


# Update from 2024 study

We have now uploaded all introners identified in our 2024 study. See file all_introners_2024.tar.gz . Introners are in exon-intron-exon fasta format. Each introner was extracted along with 20bp exoninc sequence on each flank. Flanking sequence is in capital letters and intronic sequence is in lowercase letters. For example:
if the sequence is AGTAGAGAGAAGTAGAGAGAgtnnnnnnnagAGTAGAGAGAAGTAGAGAGA, gtnnnnnnnag is the intron.

