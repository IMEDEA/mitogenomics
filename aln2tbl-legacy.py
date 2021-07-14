#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# This software is released under the license GNU GPLv3
#title           :fastalg2tbl.py
#description     :convert an assembly in fasta format into a feature table format
#author          :J. PONS, J.J. Enseñat
#date            :20190514
#version         :0.1
#usage           :fastalg2tbl.py -i infile.fasta -p file.txt -table 5 > outfile.tbl
#notes           :
#python_version  :2.7
#comments:
# USAGE DEFINITION
# "EXTREMELY IMPORTANT!!!!!!!!!"
# "HOW TO NAME GENES IN FASTA INPUT FILE"
# "NAME tRNA coding genes using single letter code e.g. A = Ala, C = Cys, EXCEPT FOR S1 = Serine1, S2 = Serine2 AND LEUCINE L1 AND L2"
# "NAME ribosomal coding genes as rrnL (or 16S) and rrnS (or 12S)"
# "NAME control region also known as A+T rich region, D-loop or hiper-variable region as CR"
# "NAME protein coding genes (CDS) as follows: atp6,atp8,cox1,cox2,cox3,cob,nad1,nad2,nad3,nad4,nad4L,nad5,nad6"
#==============================================================================


import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(prog="fastal2tbl.py",usage="fastalg2tbl.py -i infile.fasta -p file.txt -table 5 > outfile.tbl",description="Description: convert an assembly in fasta format into a feature table format")
parser.add_argument("-i", required=True, help="file with assembly in fasta format")
parser.add_argument('-p', required=True, help='Genes encoded in positive strand: file with gene names in a single line and separated by commas')
parser.add_argument('-table', type=int, required=True, help='Number of the Genetic Code Translation Table: e..g. Mito Inv table 5')
args = parser.parse_args()

#LISTS AND DICTIONARIES DEFINITIONS

#Define list of CDS genes (protein coding genes) in mitochondrial genome
c = ["atp6","atp8","cox1","cox2","cox3","cob","nad1","nad2","nad3","nad4","nad4L","nad5","nad6"]

#Define list of tRNA genes in mitochondrial genome
t = ["A","C","D","E","F","G","H","I","K","L1","L2","M","N","P","Q","R","S1","S2","T","V","W","Y"]

#Define list of ribosomal genes in mitochondrial genome
r = ["rrnL","rrnS","16S","12S"]

#Define control region (AT rich region) in mitochondrial genome
#use CR1 and CR2 in mitochondrial genomes with two control regions
n = ["CR","control_region","CR1","CR2"]

#Define starts and stop codons for different genetic codes (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
if (args.table) == 2:
	i = ["ATT","ATC","ATA","ATG","GTG"] #Vertebrate Mitochondrial Code (table=2) CANONICAL START CODONS
	s = ["TAA","TAG","AGA","AGG"] #Vertebrate Mitochondrial Code (table=2) CANONICAL STOP CODONS

elif (args.table) == 3:
	i = ["ATA","ATG","GTG"] #Yeast (table=3) Mitochondrial Code CANONICAL START CODONS
	s = ["TAA","TAG"] #Yeast (table=3) Mitochondrial Code CANONICAL STOP CODONS

elif (args.table) == 4:
	i = ["TTA","TTG","CTG","ATT","ATC","ATA","ATG","GTG"] #Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma (table=4) and Echinoderm and Flatworm Mitochondrial Code (table=9) CANONICAL START CODONS
	s = ["TAA","TAG"] #Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma (table=4) and Echinoderm and Flatworm Mitochondrial Code (table=9) CANONICAL STOP CODONS 

elif (args.table) == 5:
	i = ["TTG","ATT","ATC","ATA","ATG","GTG"] #Invertebrate (table=5) Mitochondrial Code CANONICAL START CODONS
	s = ["TAA","TAG"] #Invertebrate (table=5) Mitochondrial Code CANONICAL STOP CODONS

elif (args.table) == 9:
	i = ["ATG","GTG"] #Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma (table=4) and Echinoderm and Flatworm Mitochondrial Code (table=9) CANONICAL START CODONS
	s = ["TAA","TAG"] #Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma (table=4) and Echinoderm and Flatworm Mitochondrial Code (table=9) CANONICAL STOP CODONS 

elif (args.table) == 13:
	i = ["TTG","ATA","ATG","GTG"] #Yeast (table=3), Invertebrate (table=5), and Ascidian (table=13) Mitochondrial Code CANONICAL START CODONS
	s = ["TAA","TAG"] #Yeast (table=3), Invertebrate (table=5), and Ascidian (table=13) Mitochondrial Code CANONICAL STOP CODONS

else:
	print "ERROR: Genetic Code Table",str(args.table),"not implemented"
	sys.exit(0)

#print "Genetic Code implemented Table",str(args.table),i,s

#Import, create and define list of genes coded in the positive (or plus or ) strand
#Take command argument -p (file with gene names in a single line and separated by commas)
positive = (args.p)
for line in open(positive):
	plus = line.strip('\n\r').split(',')
#print str(plus) #Uncomment this line to print out which genes were included as being coded on the positive strand

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

#IMPORT ALL ARGUMENTS FROM TERMINAL
filename = (args.i)
fasta_sequences = SeqIO.parse(open(filename),'fasta')

#PRINT THE NAME OF THE SPECIE OR MITOGENOME TO BE ANNOTATED (GENERALLY FIRST SEQUENCE WHOLE GENOME ON WHERE FEATURES WILL BE ALIGNED/MAPPED
title = list(SeqIO.parse(open(filename),'fasta'))
print ">Features",(title[0].id)

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
			print str(pos)+"\t"+str(pos_trailing)+"\tgene"
			print "\t\t\tgene"+"\t"+str(seq_record.id)
			print str(pos)+"\t"+str(pos_trailing)+"\tCDS"
			print "\t\t\tproduct"+"\t"+str(mydic[seq_record.id])
			startbaseone = pos-1
			startbasethree = pos+2
			startingcodon = seq_record.seq[startbaseone:startbasethree]
			ataatg = str(startingcodon)
			length = pos_trailing - pos +1
#			print "length",length,"startcodon",ataatg
			if ataatg not in i: #If start codon is not canonical start codon (LIST i)
				#print "\t\t\ttransl_except"+"\t(pos:"+str(pos)+".."+str(pos+2)+","+str(ataatg)+",aa:Met)"
				print "\t\t\ttransl_except"+"\t(pos:"+str(pos)+".."+str(pos+2)+",aa:Met)"
			else: pass
			if length%3 == 1 and str(seq_record.seq[pos_trailing - 1]) == "T": #If stop codon is truncated and composed of a simple base e.g. (T)
#				print "truncated stop codon",str(seq_record.seq[position_trailing])
				print "\t\t\ttransl_except"+"\t(pos:"+str(pos_trailing)+",aa:TERM)"
				#print "\t\t\tnote"+"\t"+"TAA stop codon is completed by the addition of 3' A residues to the mRNA"
			elif length%3 == 2 and str(seq_record.seq[pos_trailing - 2:pos_trailing]) == "TA": #If stop codon truncated but composed of two bases e.g. (TA)
#				print "truncated stop codon",str(seq_record.seq[position_trailing-1:position_trailing+1])
				print "\t\t\ttransl_except"+"\t(pos:"+str(pos_trailing - 1)+".."+str(pos_trailing)+",aa:TERM)"
				#print "\t\t\tnote"+"\t"+"TAA stop codon is completed by the addition of 3' A residues to the mRNA"
			elif length%3 == 0 and str(seq_record.seq[pos_trailing-3:pos_trailing]) in s:
				pass
#				print "stop is canonical", str(seq_record.seq[position_trailing-2:position_trailing+1])
			else:
				print seq_record.id,"ERROR IN STOP CODON"
		elif seq_record.id in t: #If gene is a tRNA (LIST t) in positive strand
			print str(pos)+"\t"+str(pos_trailing)+"\tgene"
			print "\t\t\tgene"+"\t"+"trn"+str(seq_record.id)
			print str(pos)+"\t"+str(pos_trailing)+"\ttRNA"
			print "\t\t\tproduct"+"\t"+"tRNA_"+str(mydic[seq_record.id])
		elif seq_record.id in r: #If gene is a ribosomal RNA (LIST r) in positive strand
			print str(pos)+"\t"+str(pos_trailing)+"\tgene"
			print "\t\t\tgene"+"\t"+str(seq_record.id)
			print str(pos)+"\t"+str(pos_trailing)+"\trRNA"
			print "\t\t\tproduct"+"\t"+str(mydic[seq_record.id])
		elif seq_record.id in n: #if gene is a control region (LIST n) in positive strand
			print str(pos)+"\t"+str(pos_trailing)+"\t"+str(mydic[seq_record.id])
			print "\t\t\tnote"+"\t"+"putative control region"
		elif seq_record.id not in c or t or r or n:
			print str(pos)+"\t"+str(pos_trailing)+"\t"+"misc_feature"
			print "\t\t\tnote"+"\t"+str(seq_record.id)



#IF GENE IS NOT IN PLUS STRAND (e.i. not in positive strand, LIST -p) THEN MEANS THAT IS ENCODED ON THE NEGATIVE STRAND
	else:
		if seq_record.id in c: #Annotate protein coding cds in negative strand
			newseq = seq_record.seq.ungap("-") #Remove gaps
			seqrev = newseq.reverse_complement() #Make reverse complement sequence
#			print str(seq_record.seq[position:position_trailing])
			print str(pos_trailing)+"\t"+str(pos)+"\tgene"
			print "\t\t\tgene"+"\t"+str(seq_record.id)
			print str(pos_trailing)+"\t"+str(pos)+"\tCDS"
			print "\t\t\tproduct"+"\t"+str(mydic[seq_record.id])
#			print seq_record.seq
#			print newseq
#			print seqrev
#			print position_trailing+1,position+1
			rev_startingcodon = seqrev[0:3] #Find out start codon of reverse complement sequence
			rev_ataatg = str(rev_startingcodon) #Convert to string
			length = len(seqrev) #Estimate sequence length
#			print "length",length,rev_ataatg
			if rev_ataatg not in i: #If start codon is not canonical start codon (LIST i)
				#print "\t\t\ttransl_except"+"\t(pos:complement("+str(position_trailing-1)+".."+str(position_trailing+1)+"),"+str(rev_ataatg)+"),aa:Met)"
				print "\t\t\ttransl_except"+"\t(pos:complement("+str(pos_trailing-2)+".."+str(pos_trailing)+"),aa:Met)"
			else:
#				print "startcodon is canonical (pos:complement(",str(position_trailing-1)+"..",str(position_trailing+1),",",str(rev_ataatg)+"),aa:Met)"
				pass
			if length%3 == 1 and seqrev[-1:] == "T": #If stop codon is truncated and composed of a simple base (T)
#				print "truncated stop codon",seqrev[-1:]
				print "\t\t\ttransl_except"+"\t(pos:complement("+str(pos)+"),aa:TERM)"
				print "\t\t\tnote"+"\t"+"TAA stop codon is completed by the addition of 3' A residues to the mRNA"
			elif length%3 == 2 and seqrev[-2:] == "TA": #If stop codon truncated but composed of two bases (TA)
#				print "truncated stop codon",seqrev[-2:]
				print "\t\t\ttransl_except"+"\t(pos:complement("+str(pos)+".."+str(pos + 1)+"),aa:TERM)"
				print "\t\t\tnote"+"\t"+"TAA stop codon is completed by the addition of 3' A residues to the mRNA"
			elif length%3 == 0 and seqrev[-3:] in s:
#				print "pos:complement",pos,pos+2,length,"stop codon is canonical",seqrev[-3:]
				pass
			else:
				print seq_record.id,"ERROR IN STOP CODON"
		elif seq_record.id in t: #If gene is a tRNA in negative strand
			print str(pos_trailing)+"\t"+str(pos)+"\tgene"
			print "\t\t\tgene"+"\t"+"trn"+str(seq_record.id)
			print str(pos_trailing)+"\t"+str(pos)+"\ttRNA"
			print "\t\t\tproduct"+"\t"+"tRNA_"+str(mydic[seq_record.id])
		elif seq_record.id in r: #If gene is a ribosomal RNA (rRNA) in negative strand
			print str(pos_trailing)+"\t"+str(pos)+"\tgene"
			print "\t\t\tgene"+"\t"+str(seq_record.id)
			print str(pos_trailing)+"\t"+str(pos)+"\trRNA"
			print "\t\t\tproduct"+"\t"+str(mydic[seq_record.id])
		elif seq_record.id in n: #If control region is in negative strand
			print str(pos_trailing)+"\t"+str(pos)+"\t"+str(mydic[seq_record.id])
			print "\t\t\tnote"+"\t"+"putative control region"

