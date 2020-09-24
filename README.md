# Provided versions

There are two scripts 

aln2tbl-2.py  -> For use with python 2. This will not be maintained.

aln2tbl-3.py  -> For use with python 3.

This software is provided as is without warranty of any kind.

# Dependencies

This software uses Biopython

```
pip install biopython
```

# Usage

Simply copy the script to your executable path and use it as follows

```
aln2tbl-3.py -f infile.fasta -g genesfile.txt -t 5 > outfile.tbl
```

Parameters:

**-f, --fasta** -> Input file with assemby in fasta format.

**-g, --genes** -> Genes encoded in positive strand: file with gene names in a single line and separated by commas.

**-t, --table** -> Number of the Genetic Code Translation Table: e..g. Mito Inv table 5.

