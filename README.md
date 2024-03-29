
# Provided versions

There are three scripts: 

aln2tbl.py -> Converts an assembly of individual genes mapped to the mitogenome sequence in fasta format into a feature table (python 3 version).

aln2tbl-legacy.py ->  python 2 version of the previous script - (this will not be maintained).

mitos2fasta.py -> Maps individual genes to a mitogenome sequence (python 3 version).

**This software is released under the license GNU GPLv3.**

**This software is provided as is without warranty of any kind.**

# Download

Download scripts and examples from https://github.com/IMEDEA/mitogenomics

In linux terminal

```
git clone https://github.com/IMEDEA/mitogenomics
```

If git is not installed visit https://github.com/git-guides/install-git  
which is available for any OS

Make scripts executable

```
cd mitogenomics
```

```
chmod +x aln2tbl.py mitos2fasta.py
```

# Dependencies

This software uses Biopython and argparse modules in python 3 (tested version 1.78-2). They can be installed in your system using pip:

```
pip install biopython argparse
```

# mitos2fasta.py


**DESCRIPTION:**

Description: align (map) genes to mitogenome sequence to build an assembly.

Type on terminal mitos2fasta.py -h for further information.

People wishing to contribute to the software, report issues or seek support can contact Joan Pons at jpons@imedea.uib-csic.es


**USAGE:**


**mitos2fasta.py -m mitofile.fas -g genesfile.fas -c convertfile_Y/N > assembly.fas**


**PARAMETERS:**

**-m MITOFILE --mitofile** -> Input file with mitogenome sequence in fasta format, as submitted to MITOS2

**-g GENESFILE, --genesfile** -> Input file with MITOS2 output with individual genes in fasta format

**-c CONVERTFILE, --convertfile** -> Gene names (fasta headers) from MITOS2 will be simplified and made compliant with aln2tbl. Yes=Y No=N

Copy and paste the next example in your terminal (~/mitogenomics)

```
./mitos2fasta.py -m ./example/input/Hyalella_solida_mitogenome.fas -g ./example/input/Hyalella_solida_genes_mitos2.fas -c Y > ./example/input/Hyalella_solida_assembly.fas
```

# aln2tbl.py

**DESCRIPTION:**

Convert an assembly in fasta format into a feature table format.

This python script was designed for those modifying manually the genes mapped on the mitochondrial genome sequence and then
automatically generate a new table fatures with updated 5' and 3' ends. It's a perfect tool to correct annotations from MITOS2
and other similar software. The new table features can be used to build a sqn file to submit annotated mitogenome to any Data Bank (ENA, NCBI and DDBJ).

Type on terminal aln2tbl.py -h for further information.

**We encourage users to run the full pipeline example (see below) on the Hyalella data provided in order to test their setup and familiarize with the commands.**

People wishing to contribute to the software, report issues or seek support can contact Joan Pons at jpons@imedea.uib-csic.es

**USAGE:**

Simply copy the script to your executable path and use it as follows


**aln2tbl.py -f assembly_file.fas -g forward_genes_file.txt -c number_genetic_code > feature_table_file.tbl**


**PARAMETERS:**

**-f, --fasta** -> input text file with assembly in fasta format

**-g, --genes** -> input text file with list of gene names coded in forward strand (plus or positive strand) in a single line and separated by commas.

**-c, --code** -> Number (integer) of the appropriate mitochondrial Genetic Code Translation Table: vertebrate (2), yeast (3), mold, protozoan and coelenterate (4), invertebrate (5), echinoderm and flatworm (9), ascidian (13)

**PARAMETERS EXPLANATION:**

**-f** -> A fasta alignment (contig or assembly) file (--fasta) including the complete nucleotide sequence of the mitogenome
as well as the sequence of each mapped gene (generally 37), with a single line per gene sequence.

 -> See example in /example/input/Hyalella_solida_assembly_manually_curated.fas

See FAQs for the correct naming of genes.

**-g** -> A plain text file (--genes) listing genes that are encoded in forward orientation with respect to the sequence being submitted, separated by commas and without spaces.
This information allows the script to correctly reverse/complement the nucleotide sequence of genes encoded on the opposite strand prior to annotation.
Control regions, if any (see FAQs below), should be also included in the gene file (--genes), since they are generally annotated in forward orientation.

 -> See example file in /example/input/forward_genes.txt

**-c** -> The last argument (--code) parses the number of the appropriated mitochondrial genetic code
(e.g. 2 for vertebrate mitochondrial code, 5 for invertebrate mitochondrial code).

Visit this web page for a complete list of genetic codes https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
See FAQs below for a summary of start and stop codons of different mitochondrial genetic codes.

**The example script can be run using the following example:**

Copy and paste the next example in your terminal (~/mitogenomics)

```
./aln2tbl.py -f ./example/input/Hyalella_solida_assembly_manually_curated.fas -g ./example/input/forward_genes.txt -c 5 > ./example/input/Hyalella_solida_feature_table.tbl
```

The output file is saved as plain text with tbl extension to be readily identified as a feature table. 

Compare your Hyalella_solida_feature_table.tbl with our output file /example/output/Hyalella_solida_feature_table.tbl


# Full pipeline example (example files in /example/)

**Input files are provided in the ~/mitogenomics/example/input folder**

Hyalella_solida_mitogenome.fas -> Mitogenome sequence in fasta format

Hyalella_solida_genes_mitos2.fas -> MITOS2 output, automatic annotation of the genome.

Hyalella_solida_assembly_manually_curated.fas -> MITOS2 output manually curated in an alignment editor (SeaView or Aliview). In this example, duplicated genes have been removed or merged, the boundaries of some genes have been corrected and names of gene names have been made compliant (see FAQs)

forward_genes.txt -> List of genes in forward orientation

submission_template.sbt -> Submission information created using the NCBI web function at https://submit.ncbi.nlm.nih.gov/genbank/template/submission/

**Pipeline steps**

1) Send mitogenome sequence "/example/input/Hyalella_solida_mitogenome.fas" to MITOS2
(http://mitos2.bioinf.uni-leipzig.de/index.py)

Select correct genetic code (Invertebrate for this example) and in advanced setting unselect OH, OL, and intron features.
Also unselect circular since this example is not a complete mitogenome as we failed to sequence part of control region. 

2) From MITOS2 output download fasta file sequences of annotated genes (FAS file link). Rename file to "Hyalella_solida_genes_mitos2.fas"

3) Maps the automatically annotated genes (e.g. file .fas from MITOS2 automatic annotation; see http://mitos2.bioinf.uni-leipzig.de/index.py) to the mitogenome to produce a fasta assembly with mitos2fasta python script

Copy and paste the next example in your terminal (~/mitogenomics)

```
./mitos2fasta.py -m ./example/input/Hyalella_solida_mitogenome.fas -g ./example/input/Hyalella_solida_genes_mitos2.fas -c Y > ./example/input/Hyalella_solida_assembly.fas
```

This same file/step can be produced using alternative strategies (e.g. bowtie) or manually.

4) If necessary (generally it is so) manually check and modify 5' and 3' ends of genes and annotate and add control regions in Seaview or AliView as we did manually. Our curated files is in:

"/example/input/Hyalella_solida_assembly_manually_curated.fas"

5) Creates feature table file with aln2tbl python script from the manually curated fasta assembly.

Copy and paste the next example in your terminal (~/mitogenomics)

```
./aln2tbl.py -f ./example/input/Hyalella_solida_assembly_manually_curated.fas -g ./example/input/forward_genes.txt -c 5 > ./example/input/Hyalella_solida_feature_table.tbl
```

6) build .sqn file to submit annotated mitogenome to GenBank/ENA and check/validade for errors

Copy and paste the next example in your terminal (~/mitogenomics)

```
tbl2asn -i ./example/input/Hyalella_solida_mitogenome.fas -f ./example/input/Hyalella_solida_feature_table.tbl -t ./example/input/submission_template.sbt -a s -V bv -T -j "[mgcode=5] [location=mitochondrion] [organism=Hyalella solida]"
```

Submission files are created using the NCBI tbl2asn script. This last commands assumes the tbl2asn script is correctly installed in your system. Please check the NCBI website for additional infromation on tbl2asn https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/

Source modifiers (organism, genetic code, ...), here provided to tbl2asn using the -j option, can be alternatively specified in the header of the genome fasta. For further information about source qualifiers (-j) see: https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/


**Output files comparison**

After running full pipeline commands, 6 additional output files will be found in ~/mitogenomics/example/input folder. 

These same files are provided in folder ~/mitogenomics/example/output for comparison.

Hyalella_solida_assembly.fas -> Genes mapped over the genome, in fasta format. This file can be edited using an alignment editor for manual curation of annotations. Sample manual curation is here provided as input file Hyalella_solida_assembly_manually_curated.fas. 
 
Hyalella_solida_feature_table.tbl -> Feature table with gene annotations.        

Hyalella_solida_mitogenome.val -> Validation output, should contain no (major) error.

errorsummary.val -> Validation output, should contain no (major) error.

Hyalella_solida_mitogenome.sqn -> Annotation of the mitogenome, ready for submission.

Hyalella_solida_mitogenome.gbf -> Human readable version of the submission file, in GenBank format.

# FAQs/FAQs/FAQs/:

**1) The script is apparently missing the correct python interpreter**

This script assumes the python3 interpreter to be in /usr/bin/env python3. If this is not the case, or multiple python
installations are available, the full path to python3 interpreter can be added to the command line
(e.g. ./usr/bin/python3 aln2tbl.py or ~/miniconda3/bin/python ./aln2tbl.py) or path exported in .bashrc file.

**2) How should individual genes be named in the input file?**

Gene names must comply with the names proposed by Boore and Brown (2000) as in most recent GenBank annotations.

Protein coding gene names (CDS) are as follows: atp6, atp8, cob, cox1, cox2, cox3, nad1, nad2, nad3, nad4, nad4L, nad5 and nad6.  

Accepted names for ribosomal genes are rrnL (or 16S) and rrnS (or 12S) for the large and small ribosomal subunit, respectively.  

tRNAs are indicated using the single letter corresponding to the encoded tRNA (e.g. M for methionine).  

The two tRNA genes that are generally present for Leucine and Serine must be differentiated by post-pending the correct
number to the gene name following the same convention used in MITOS: L1 for tRNALeu(CUN), L2 for tRNALeu(UUR), S1 for tRNASer(AGN) and S2 for tRNASer(UCN)). 

The control region, if included, must be named 'control_region' or 'CR'. Up to two codon regions can be annotated ‘CR1’ and ‘CR2’.

**3) Should all spacers be annotated as Control Regions?**

No, do not annotate any non-coding region as a Control region. The Control region, just like PCGs, rRNAs and tRNAs, should be identified based on positive evidence (position, base composition, presence of specific conserved elements) prior to annotation. 

**4) Which start and stop codons are recognized as standard in each genetic code?**

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


**5) How should I obtain the starting file (i.e. an assembly of individual genes mapped to the mitogenome sequence)?**

The initial fasta alignment of automatically annotated gene sequences against the complete mitogenome can be obtained
in several ways starting from the results of an automatic annotation (e.g. MITOS2). One flexible option, using free software
that can be called using bash scripting, is to build the contig using Bowtie2 (Langmead & Salzberg, 2012), sort the bam file
by coordinate in picard tools (http://broadinstitute.github.io/picard), and finally convert the bam file to fasta with the
python3 script sam2fasta.py (https://sourceforge.net/projects/sam2fasta/files/).

Another option, provided that gene and genome sequences are position ordered, identical and ungapped (e.g. in the MITOS2 output),
is to use the python 3 script mitos2aln.py, provided as companion to aln2tbl.phy. This script produces a fasta alignment of
each gene nucleotide sequence relative to complete mitogenome sequence and parse mitos gene nomenclature of fasta headers
with four semicolon separated fields to single name compatible with aln2tbl.py naming. 

In both options, the fasta contig can be the easily refined manually in a visual alignment editor such as SeaView or Aliview
prior to further processing in aln2tbl.py. 

**6) Is it possible to annotated a TA truncated stop codon?**

Both T and TA truncated stop codons are correctly processed by aln2tbl.py.
Nevertheless the use of truncated TA stop codons is discouraged as this is not recognized by downstream application tbl2asn.
A simple workaround is to indicate a T truncated stop codon in any case.

**7) Can the script handle partial genes?**

No. The aln2tbl.py script can not handle partial genes, though this is an uncommon circumstance in mitogenomes.
If present, the feature table must be manually corrected by adding the symbols < and > for partial 5’ and 3’ ends, respectively.

See link for futher information 
https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/#partial_CDS

**8) Can the script handle duplicated tRNAs (apart from Leucine and Serine)?**

No. The presence of duplicated tRNA genes, which is also uncommon, cannot be handled as it results in a duplicated gene name. 
These must be manually annotated in the feature table.

**9) Can the script handle genes including first position as internal one e.g. cox1 14399-1420?**

No. If the genome is circular, we foster the good practice of linearizing the sequence at the beginning of one gene, generally tRNA-I to circumvent the problem and produce a well organized submission. 

**10) The genome is incomplete and missing an inner portion of a single gene**

This occurrence must be manually annotated in the feature table.

