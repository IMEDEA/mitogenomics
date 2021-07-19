#!/usr/bin/env python3
# This software is released under the license GNU GPLv3

import sys
import argparse
from Bio import SeqIO

# USAGE DEFINITION
# "EXTREMELY IMPORTANT!!!!!!!!!"
# "HOW TO NAME GENES IN FASTA INPUT FILE"
# "NAME tRNA coding genes using single letter code e.g. A = Ala, C = Cys, EXCEPT FOR S1 = Serine1, S2 = Serine2 AND LEUCINE L1 AND L2"
# "NAME ribosomal coding genes as rrnL (or 16S) and rrnS (or 12S)"
# "NAME control region also known as A+T rich region, D-loop or hiper-variable region as CR"
# "NAME protein coding genes (CDS) as follows: atp6,atp8,cox1,cox2,cox3,cob,nad1,nad2,nad3,nad4,nad4L,nad5,nad6"

#LISTS AND DICTIONARIES DEFINITIONS

#Define list of CDS genes (protein coding genes) in mitochondrial genome
c = ["atp6", "atp8", "cox1", "cox2", "cox3", "cob", "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]

#Define list of tRNA genes in mitochondrial genome
t = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L1", "L2", "M", "N", "P", "Q", "R", "S1", "S2", "T", "V", "W", "Y"]

#Define list of ribosomal genes in mitochondrial genome
r = ["rrnL", "rrnS", "16S", "12S"]

#Define control region (AT rich region) in mitochondrial genome
#use CR1 and CR2 in mitochondrial genomes with two control regions. You can add up to 10 control regions
n = ["CR", "control_region", "CR1", "CR2"]

#DEFINE DICTIONARY CORRESPONDENCE BETWEEN "GENE" : "PRODUCT"
mydic = {
"A" : "Ala", 
"C" : "Cys", 
"D" : "Asp", 
"E" : "Glu", 
"F" : "Phe", 
"G" : "Gly", 
"H" : "His", 
"I" : "Ile", 
"K" : "Lys", 
"L1" : "Leu1", 
"L2" : "Leu2", 
"M" : "Met", 
"N" : "Asn", 
"P" : "Pro", 
"Q" : "Gln", 
"R" : "Arg", 
"S1" : "Ser1", 
"S2" : "Ser2", 
"T" : "Thr", 
"V" : "Val", 
"W" : "Trp", 
"Y" : "Tyr", 
"atp6" : "ATP synthase F0 subunit 6", 
"atp8" : "ATP synthase F0 subunit 8", 
"cox1" : "cytochrome c oxidase subunit 1", 
"cox2" : "cytochrome c oxidase subunit 2", 
"cox3" : "cytochrome c oxidase subunit 3", 
"cob" : "cytochrome b", 
"nad1" : "NADH dehydrogenase subunit 1", 
"nad2" : "NADH dehydrogenase subunit 2", 
"nad3" : "NADH dehydrogenase subunit 3", 
"nad4" : "NADH dehydrogenase subunit 4", 
"nad4L" : "NADH dehydrogenase subunit 4L", 
"nad5" : "NADH dehydrogenase subunit 5", 
"nad6" : "NADH dehydrogenase subunit 6", 
"rrnL" : "16S ribosomal RNA", 
"16S" : "16S ribosomal RNA", 
"rrnS" : "12S ribosomal RNA", 
"12S" : "12S ribosomal RNA", 
"CR" : "misc_feature", 
"CR1" : "misc_feature", 
"CR2" : "misc_feature", 
"control_region" : "misc_feature"
}


def get_args():
    parser = argparse.ArgumentParser(prog = "aln2tbl.py", usage = "aln2tbl.py -f assembly_file.fas -g forward_genes.txt -c number_genetic_code > feature_table_file.tbl", description = "Description: convert an assembly in fasta format into a feature table format. Type on terminal aln2tbl.py -h for further information. For additional information visit https://github.com/IMEDEA/mitogenomics or email Joan Pons at jpons@imedea.uib-csic.es")
    parser.add_argument("-f", "--fasta", required = True, help = "input text file with assembly in fasta format: example file in /example/output/Hyalella_solida2319A_assembly_manually_curated.fas")
    parser.add_argument("-g", "--genes", required = True, help = "input text file with list of gene names coded in forward strand (plus or positive strand) in a single line and separated by commas: example file in /example/input/forward_genes.txt")
    parser.add_argument("-c", "--code", type = int, required = True, help = "Number (integer) of the appropriate mitochondrial Genetic Code Translation Table: vertebrate (2), yeast (3), mold, protozoan and coelenterate (4), invertebrate (5), echinoderm and flatworm (9), ascidian (13)")
    args = parser.parse_args()
    
    #Define starts and stop codons for different genetic codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
    if args.code == 2:
        i = ["ATT","ATC","ATA","ATG","GTG"] #Vertebrate Mitochondrial Code (genetic_code=2) CANONICAL START CODONS
        s = ["TAA","TAG","AGA","AGG"] #Vertebrate Mitochondrial Code (genetic_code=2) CANONICAL STOP CODONS

    elif args.code == 3:
        i = ["ATA","ATG","GTG"] #Yeast (genetic_code=3) Mitochondrial Code CANONICAL START CODONS
        s = ["TAA","TAG"] #Yeast (genetic_code=3) Mitochondrial Code CANONICAL STOP CODONS

    elif args.code == 4:
        i = ["TTA","TTG","CTG","ATT","ATC","ATA","ATG","GTG"] #Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma (table=4) CANONICAL START CODONS
        s = ["TAA","TAG"] #Mold, Protozoan, and Coelenterate Mitochondrial Code (genetic_code=4) CANONICAL STOP CODONS 

    elif args.code == 5:
        i = ["TTG","ATT","ATC","ATA","ATG","GTG"] #Invertebrate (table=5) Mitochondrial Code CANONICAL START CODONS
        s = ["TAA","TAG"] #Invertebrate (table=5) Mitochondrial Code CANONICAL STOP CODONS

    elif args.code == 9:
        i = ["ATG","GTG"] #Echinoderm and Flatworm Mitochondrial Code (genetic_code=9) CANONICAL START CODONS
        s = ["TAA","TAG"] #Echinoderm and Flatworm Mitochondrial Code (genetic_code=9) CANONICAL STOP CODONS 

    elif args.code == 13:
        i = ["TTG","ATA","ATG","GTG"] #Ascidian (genetic_code=13) Mitochondrial Code CANONICAL START CODONS
        s = ["TAA","TAG"] #Ascidian (genetic_code=13) Mitochondrial Code CANONICAL STOP CODONS
    else:
        print("ERROR: Genetic Code Table {} not implemented".format(str(args.table)))
        sys.exit(-1)

    #Import, create and define list of genes coded in the positive (or plus or ) strand
    #Take command argument -p (file with gene names in a single line and separated by commas)
    genesfile = args.genes
    try:
        with open(genesfile) as gf:
            plus = gf.readline().strip('\n\r').split(',')
    except:
        print("ERROR: Genes file cannot be read")
        sys.exit(-2)


    #IMPORT ALL ARGUMENTS FROM TERMINAL
    fastafile = args.fasta
    try:
        fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    except:
        print("ERROR: Fasta file cannot be read")
        sys.exit(-3)
    
    #PRINT THE NAME OF THE SPECIE OR MITOGENOME TO BE ANNOTATED (GENERALLY FIRST SEQUENCE WHOLE GENOME ON WHERE FEATURES WILL BE ALIGNED/MAPPED
    title = list(SeqIO.parse(open(fastafile),'fasta'))
    print(">Features {}".format(title[0].id))
    
    return i, s, plus, fasta_sequences


def process(i, s, plus, fasta_sequences):
    
    #ESTIMATE 5' AND 3' END OF EACH GENE FEATURE FROM FASTA ASSEMBLY/ALIGNMENT/CONTIG/MAPPING
    #STRATEGY: ESTIMATE POSITION BY COUNTING LEADING AND TRAILING GAPS SINCE FEATURES ARE ALIGNED/MAPPED TO THE WHOLE MITOCHONDRIAL GENOME
    
    for seq_record in fasta_sequences:
        seq_record.seq = seq_record.seq.upper() #this makes sequences upper case using function in Bio import SeqIO

        ltrimmed = seq_record.seq.lstrip('-')
        ttrimmed = seq_record.seq.strip('-')
        pos = len(seq_record) - len(ltrimmed) + 1
        pos_trailing = pos + len(ttrimmed) - 1

        #START GENE ANNOTATION IN TABLE FORMAT TO THEM BE READED IN tbl2asn AND GET GENBANK OR EMBL FORMAT TO SUBMIT THEM
        if seq_record.id in plus: #If gene is coded on positive strand (LIST -p)
            if seq_record.id in c: #If gene is a protein coding cds (LIST c)
                print("{}\t{}\tgene".format(str(pos),str(pos_trailing)))
                print("\t\t\tgene\t{}".format(str(seq_record.id)))
                print("{}\t{}\tCDS".format(str(pos),str(pos_trailing)))
                print("\t\t\tproduct\t{}".format(str(mydic[seq_record.id])))
                startbaseone = pos - 1
                startbasethree = pos + 2
                startingcodon = seq_record.seq[startbaseone:startbasethree]
                ataatg = str(startingcodon)
                length = pos_trailing - pos + 1
                ##print "length",length,"startcodon",ataatg
                if ataatg not in i: #If start codon is not canonical start codon (LIST i)
                    print("\t\t\ttransl_except\t(pos:{}..{},aa:Met)".format(str(pos), str(pos + 2)))
                else: 
                    pass
                if length%3 == 1 and str(seq_record.seq[pos_trailing - 1]) == "T": #If stop codon is truncated and composed of a simple base e.g. (T)
                    print("\t\t\ttransl_except\t(pos:{},aa:TERM)".format(str(pos_trailing)))
                    #print("\t\t\tnote\tTAA stop codon is completed by the addition of 3' A residues to the mRNA")
                elif length%3 == 2 and str(seq_record.seq[pos_trailing - 2:pos_trailing]) == "TA": #If stop codon truncated but composed of two bases e.g. (TA)
                    print("\t\t\ttransl_except\t(pos:{}..{},aa:TERM)".format(str(pos_trailing - 1), str(pos_trailing)))
                    #print("\t\t\tnote\tTAA stop codon is completed by the addition of 3' A residues to the mRNA")
                elif length%3 == 0 and str(seq_record.seq[pos_trailing-3:pos_trailing]) in s:
                    pass
                else:
                    print("{} ERROR IN STOP CODON".format(seq_record.id))
            elif seq_record.id in t: #If gene is a tRNA (LIST t) in positive strand
                print("{}\t{}\tgene".format(str(pos), str(pos_trailing)))
                print("\t\t\tgene\ttrn{}".format((seq_record.id)))
                print("{}\t{}\ttRNA".format(str(pos), str(pos_trailing)))
                print("\t\t\tproduct\ttRNA_{}".format(str(mydic[seq_record.id])))
            elif seq_record.id in r: #If gene is a ribosomal RNA (LIST r) in positive strand
                print("{}\t{}\tgene".format(str(pos), str(pos_trailing)))
                print("\t\t\tgene\t{}".format(str(seq_record.id)))
                print("{}\t{}\trRNA".format(str(pos), str(pos_trailing)))
                print("\t\t\tproduct\t{}".format(str(mydic[seq_record.id])))
            elif seq_record.id in n: #if gene is a control region (LIST n) in positive strand
                print(str(pos)+"\t"+str(pos_trailing)+"\t"+str(mydic[seq_record.id]))
                print("\t\t\tnote\tputative control region")
            elif seq_record.id not in c or t or r or n:
                print("{}\t{}\tmisc_feature".format(str(pos), str(pos_trailing)))
                print("\t\t\tnote\t{}".format(str(seq_record.id)))

        #IF GENE IS NOT IN PLUS STRAND (e.i. not in positive strand, LIST -p) THEN MEANS THAT IS ENCODED ON THE NEGATIVE STRAND
        else:
            if seq_record.id in c: #Annotate protein coding cds in negative strand
                newseq = seq_record.seq.ungap("-") #Remove gaps
                seqrev = newseq.reverse_complement() #Make reverse complement sequence
                print("{}\t{}\tgene".format(str(pos_trailing), str(pos)))
                print("\t\t\tgene\t{}".format(str(seq_record.id)))
                print("{}\t{}\tCDS".format(str(pos_trailing), str(pos)))
                print("\t\t\tproduct\t{}".format(str(mydic[seq_record.id])))
                rev_startingcodon = seqrev[0:3] #Find out start codon of reverse complement sequence
                rev_ataatg = str(rev_startingcodon) #Convert to string
                length = len(seqrev) #Estimate sequence length
                if rev_ataatg not in i: #If start codon is not canonical start codon (LIST i)
                    print("\t\t\ttransl_except\t(pos:complement({}..{}),aa:Met)".format(str(pos_trailing - 2), str(pos_trailing)))
                else:
                    pass
                if length%3 == 1 and seqrev[-1:] == "T": #If stop codon is truncated and composed of a simple base (T)
                    print("\t\t\ttransl_except\t(pos:complement({}),aa:TERM)".format(str(pos)))
                    print("\t\t\tnote\tTAA stop codon is completed by the addition of 3' A residues to the mRNA")
                elif length%3 == 2 and seqrev[-2:] == "TA": #If stop codon truncated but composed of two bases (TA)
                    print("\t\t\ttransl_except\t(pos:complement({}..{}),aa:TERM)".format(str(pos), str(pos + 1)))
                    print("\t\t\tnote\tTAA stop codon is completed by the addition of 3' A residues to the mRNA")
                elif length%3 == 0 and seqrev[-3:] in s:
                    pass
                else:
                    print("{} ERROR IN STOP CODON".format(seq_record.id))
            elif seq_record.id in t: #If gene is a tRNA in negative strand
                print("{}\t{}\tgene".format(str(pos_trailing), str(pos)))
                print("\t\t\tgene\ttrn{}".format(str(seq_record.id)))
                print("{}\t{}\ttRNA".format(str(pos_trailing), str(pos)))
                print("\t\t\tproduct\ttRNA_{}".format(str(mydic[seq_record.id])))
            elif seq_record.id in r: #If gene is a ribosomal RNA (rRNA) in negative strand
                print("{}\t{}\tgene".format(str(pos_trailing), str(pos)))
                print("\t\t\tgene\t{}".format(str(seq_record.id)))
                print("{}\t{}\trRNA".format(str(pos_trailing), str(pos)))
                print("\t\t\tproduct\t{}".format(str(mydic[seq_record.id])))
            elif seq_record.id in n: #If control region is in negative strand
                print("{}\t{}\t{}".format(str(pos_trailing), str(pos), str(mydic[seq_record.id])))
                print("\t\t\tnote\tputative control region")


def main():
    i, s, plus, fasta_sequences = get_args()
    
    process(i, s, plus, fasta_sequences)


if __name__ == "__main__":
    main()
