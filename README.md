
# Provided versions

There are three scripts 

aln2tbl.py -> Converts an assembly of individual genes mapped to the mitogenome sequence in fasta format into a feature table (python 3 version).

aln2tbl-legacy.py  -> Maps individual genes to a mitogenome sequence (pyhton 2 - this will not be maintained).

# This software is released under the license GNU GPLv3

# This software is provided as is without warranty of any kind.

# Dependencies

This software uses Biopython (tested version 1.78-2)

```
pip install biopython
```

# aln2tbl

This python script was designed for those modifying manually the genes mapped on the mitochondrial genome sequence and then
automatically generate a new table fatures with updated 5' and 3' ends. It's a perfect tool to correct annotations from MITOS2
and other similar software.

People wishing to contribute to the software, report issues or seek support can contact Joan Pons at jpons@imedea.uib-csic.es

#Usage

Simply copy the script to your executable path and use it as follows

```
aln2tbl.py --fasta My_assembly.fas --genes My_forward_genes.txt --code integer > My_feature_table.tbl
```

Parameters:

**-f, --fasta** -> input text file with assembly in fasta format: example file in /mitogenomics-master/data_examples/Hyalella_solida2319A_assembly.fas

**-g, --genes** -> input text file with list of gene names coded in forward strand (plus or positive strand) in a single line and separated by commas: example file in /mitogenomics-master/data_examples/forward_genes.txt

**-c, --code** -> Number (integer) of the appropriate mitochondrial Genetic Code Translation Table: vertebrate (2), yeast (3), mold, protozoan and coelenterate (4), invertebrate (5), echinoderm and flatworm (9), ascidian (13)

The aln2table.py script takes as input

a) A fasta alignment (contig or assembly) file (--fasta) including the complete nucleotide sequence of the mitogenome
as well as the sequence of each mapped gene (generally 37), with a single line per gene sequence.

See example in /mitogenomics-master/data_examples/Hyalella_solida2319A_assembly.fas
See FAQs for the correct naming of genes.

b) A plain text file (--genes) listing genes that are encoded in forward orientation with respect to the sequence being submitted, separated by commas and without spaces.
This information allows the script to correctly reverse/complement the nucleotide sequence of genes encoded on the opposite strand prior to annotation. Control regions, if any (see FAQs below), should be also included in the gene file (--genes), since they are generally annotated in forward orientation. 

	See example in /mitogenomics-master/data_examples/forward_genes.txt

c) The last argument (--code) parses the number of the appropriated mitochondrial genetic code
(e.g. 2 for vertebrate mitochondrial code, 5 for invertebrate mitochondrial code)

Visit this web page for a complete list of genetic codes https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
See FAQs below for a summary of start and stop codons of different mitochondrial genetic codes.

The example script can be run using the following example:

```
aln2tbl.py --fasta Hyalella_solida2319A_assembly.fas --genes forward_genes.txt --code 5 > Hyalella_solida2319A.tbl
```
The output file is saved as plain text with tbl extension to be readily identified as a feature table. 


# FAQs/FAQs/FAQs/:

1) The script is apparently missing the correct python interpreter!

This script assumes the python3 interpreter to be in /usr/bin/env python3. If this is not the case, or multiple python
installations are available, the full path to python3 interpreter can be added to the command line
(e.g. ./usr/bin/python3 aln2tbl.py) or path exported in .bashrc file.

2) How should individual genes be named in the input file?
Gene names must comply with the names proposed by Boore and Brown (2000) as in most recent GenBank annotations.

Protein coding gene names (CDS) are as follows: atp6, atp8, cob, cox1, cox2, cox3, nad1, nad2, nad3, nad4, nad4L, nad5 and nad6.
Accepted names for ribosomal genes are rrnL (or 16S) and rrnS (or 12S) for the large and small ribosomal subunit, respectively.
tRNAs are indicated using the single letter corresponding to the encoded tRNA (e.g. M for methionine).
The two tRNA genes that are generally present for Leucine and Serine must be differentiated by post-pending the correct
number to the gene name following the same convention used in MITOS: L1 for tRNALeu(CUN), L2 for tRNALeu(UUR), S1 for tRNASer(AGN) and S2 for tRNASer(UCN)).
The control region, if included, must be named 'control_region' or 'CR'. Up to two codon regions can be annotated ‘CR1’ and ‘CR2’

3) Should all spacers be annotated as Control Regions?
No, do not annotate any non-coding region as a Control region. The Control region, just like PCGs, rRNAs and tRNAs, should be identified based on positive evidence (position, base composition, presence of specific conserved elements) prior to annotation. 

4) Which start and stop codons are recognized as standard in each genetic code?

According to NCBI standards, start and stop codons for different genetic codes for mitochondrial genomes are:

   Vertebrate Mitochondrial Code (genetic_code=2) :
       START CODONS for Met = ["ATT","ATC","ATA","ATG","GTG"]
       STOP CODONS = ["TAA","TAG","AGA","AGG"]

   Yeast Mitochondrial Code (genetic_code=3) :
       START CODONS for Met = ["ATA","ATG","GTG"]
       STOP CODONS = ["TAA","TAG"]

    Mold, Protozoan, and Coelenterate Mitochondrial Code (genetic_code=4) :
       START CODONS for Met = ["TTA","TTG","CTG","ATT","ATC","ATA","ATG","GTG"]
       STOP CODONS  = ["TAA","TAG"]

   Invertebrate Mitochondrial Code (genetic_code=5) :
       START CODONS for Met = ["TTG","ATT","ATC","ATA","ATG","GTG"]
       STOP CODONS = ["TAA","TAG"]

   Echinoderm and Flatworm Mitochondrial Code (genetic_code=9) :
       START CODONS for Met = ["ATG","GTG"]
       STOP CODONS  = ["TAA","TAG"]

   Ascidian Mitochondrial Code (genetic_code=13) :
       START CODONS for Met = ["TTG","ATA","ATG","GTG"]
       STOP CODONS = ["TAA","TAG"]


5) How so I obtain the starting file (i.e. an assembly of individual genes mapped to the mitogenome sequence)?
The initial fasta alignment of automatically annotated gene sequences against the complete mitogenome can be obtained
in several ways starting from the results of an automatic annotation (e.g. MITOS2). One flexible option, using free software
that can be called using bash scripting, is to build the contig using Bowtie2 (Langmead & Salzberg, 2012), sort the bam file
by coordinate in picard tools (http://broadinstitute.github.io/picard), and finally convert the bam file to fasta with the
python3 script sam2fasta.py (https://sourceforge.net/projects/sam2fasta/files/).

Another option, provided that gene and genome sequences are position ordered, identical and ungapped (e.g. in the MITOS2 output),
is to use the python v3 script mitos2aln.py, provided as companion to aln2tbl.phy. This script produces a fasta alignment of
each gene nucleotide sequence relative to complete mitogenome sequence and parse mitos gene nomenclature of fasta headers
with four semicolon separated fields to single name compatible with aln2tbl.py naming. 

In both options, the fasta contig can be the easily refined manually in a visual alignment editor such as SeaView or Aliview
prior to further processing in aln2tbl.phy. 

6) Is it possible to annotated a TA truncated stop codon?
Both T and TA truncated stop codons are correctly processed by aln2tbl.py.
Nevertheless the use of truncated TA stop codons is discouraged as this is not recognized by downstream application tbl2asn.
A simple workaround is to indicate a T truncated stop codon in any case.

7) Can the script handle partial genes?
No. The aln2tbl.py script can not handle partial genes, though this is an uncommon circumstance in mitogenomes.
If present, the feature table must be manually corrected by adding the symbols < and > for partial 5’ and 3’ ends, respectively.

8) Can the script handle duplicated tRNAs (apart from Leucine and Serine)?
No. The presence of duplicated tRNA genes, which is also uncommon, cannot be handled as it results in a duplicated gene name. 
These must be manually annotated in the feature table.

9) Can the script handle genes including first position as internal one e.g. cox1 14399-1420?
No. If the genome is circular, we foster the good practice of linearizing the sequence at the beginning of one gene, generally tRNA-I to circumvent the problem and produce a well organized submission. 

10) If the genome is incomplete and missing an inner portion of a single gene, this occurrence must be manually annotated in the feature table.



# Usage mitos2fasta

```
mitos2fasta.py -m My_mitogenome_sequence.fas -g My_mitos2_genes.fas -c Y > My_assembly.fas
```

Parameters:

**-m MITOFILE --mitofile** -> Input file with Mitogenome sequence in fasta format, as submitted to MITOS2

**-g GENESFILE, --genesfile** -> Input file with MITOS2 output with individual genes in fasta format

**-c CONVERTFILE, --convertfile** -> Gene names (fasta headers) from MITOS2 will be simplified and made compliant with aln2tbl. Yes=Y No=N

Example:
mitos2fasta.py -m Hyalella_solida2319A_mitogenome_reference_mitos2_input.fas -g Hyalella_solida2319A_genes_manually_curated.fas -c Y > Hyalella_solida2319A_assembly.fas

# Full pipeline example (example files in /mitogenomics-master/data_examples/)

1) Send mitogenome sequence "/example/input/Hyalella_solida2319A_mitogenome_reference_mitos2_input.fasta" to MITOS2
(http://mitos2.bioinf.uni-leipzig.de/index.py)

	Select correct genetic code (Invertebrate for this example) and in advanced setting unselect OH, OL, and intron features.
	Also unselect circular since this example is not a complete mitogenome since we failed to sequence part of control region. 

2) From MITOS2 output download fasta file sequences of annotated genes (FAS file link). Rename file to "./example/output/Hyalella_solida2319A_genes_mitos2_output.fas"

3) build assembly (contig) aligning/mapping genes relative to mitogenome sequence with mitos2fasta python script
```
./mitos2fasta.py -m ./example/input/Hyalella_solida2319A_mitogenome_reference_mitos2_input.fas -g ./example/output/Hyalella_solida2319A_genes_mitos2_output.fas -c Y > ./example/output/Hyalella_solida2319A_assembly.fas
```

4) If necessary (generally it is so) manually check and modify 5' and 3' ends of genes and add control regions in Seaview or AliView.
"./example/output/Hyalella_solida2319A_assembly_manually_curated.fas"

5) build feature table file with aln2tbl python script

```
./aln2tbl.py -f ./example/output/Hyalella_solida2319A_assembly_manually_curated.fas -g ./example/input/forward_genes.txt -c 5 > ./example/output/Hyalella_solida2319A.tbl
```

6) build .sqn file to submit annotated mitogenome to GenBank/ENA and check/validade for errors
```
tbl2asn -i ./example/input/Hyalella_solida2319A_linear_tbl2asn_input.fas -f ./example/output/Hyalella_solida2319A.tbl -t ./example/input/template.sbt -a s -V vb -T
```
	output files:
				file in GenBank format "./example/output/Hyalella_solida2319A_linear_tbl2asn_input.gbf"
				file in sqn format to submit to GenBank (or others repositories such as ENA and DDBJ)  "./example/output/Hyalella_solida2319A_linear_tbl2asn_input.sqn"
				file with ERROR and WARNING "./example/output/Hyalella_solida2319A_linear_tbl2asn_input.val"
				file with ERRORS "./example/output/errorsummary.val"

