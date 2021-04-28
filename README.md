# Provided versions

There are three scripts 

aln2tbl-legacy.py  -> Convert an assembly mapped genes to mitogenome sequence in fasta format into a feature table format (python 2 version - this will not be maintained).

aln2tbl.py  -> Convert an assembly mapped genes to mitogenome sequence in fasta format into a feature table format (python 3 version).

mitos2fasta.py  -> Map genes to mitogenome sequence reference (pyhton 3).


This software is provided as is without warranty of any kind.

# Dependencies

This software uses Biopython

```
pip install biopython
```

# Usage aln2tbl

Simply copy the script to your executable path and use it as follows

```
aln2tbl.py -f infile.fasta -g genesfile.txt -t 5 > outfile.tbl
```

Parameters:

**-f, --fasta** -> Input file with assemby in fasta format.

**-g, --genes** -> Genes encoded in positive strand: file with gene names in a single line and separated by commas.

**-t, --table** -> Number of the Genetic Code Translation Table: e..g. Mito Inv table 5.


# Usage mitos2fasta

```
mitos2fasta.py --mitofile reference.fasta --genefile genesfile.fasta --convertfile Y/N > alignment.fasta
```

Parameters:

**-m MITOFILE --mitofile** -> Input file with mitogenome sequence in fasta format

**-g GENESFILE, --genesfile** -> Input file with annotated genes names in fasta format

**-c CONVERTFILE, --convertfile** -> Headers of Input fasta file from mitos will be sanitaized Yes=Y No=N



