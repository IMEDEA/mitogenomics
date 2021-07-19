#!/usr/bin/env python3
# This software is released under the license GNU GPLv3

import sys
import os
import tempfile
import argparse
from Bio import SeqIO
import re


def get_args():
    parser = argparse.ArgumentParser(prog = "mitos2fasta.py", usage = "mitos2fasta.py -m mitofile.fas -g genesfile.fas -c covertfile_Y/N > assembly.fas", description = "Description: align (map) genes to mitogenome sequence. Type on terminal mitos2fasta.py -h for further information or visit https://github.com/IMEDEA/mitogenomics")
    parser.add_argument("-m", "--mitofile", required = True, help = "Input file with mitogenome sequence in fasta format, as submitted to MITOS")
    parser.add_argument("-g", "--genesfile", required = True, help = "Input file with MITOS output with individual genes in fasta format")
    parser.add_argument("-c", "--convertfile", required = True, help = "Gene names (fasta headers) from MITOS will be simplified and made compliant with aln2tbl. Yes=Y No=N")
    args = parser.parse_args()

    return args

# MAP GENES OVER MITOGENOMES SEQUENCE
def process(mitogenome, genes):
    for ref in SeqIO.parse(mitogenome, "fasta"):
        mitogenome_seq = str(ref.seq)
        # PRINT MITOGENOME
        print(">" + ref.id + "\n" + mitogenome_seq)

    for record in SeqIO.parse(genes, "fasta"):
        record_seq = str(record.seq)
        match = re.search(record_seq, mitogenome_seq, re.IGNORECASE)

        if match is None:
            # IF THERE IS NO MATCH GENE SEQUENCE IS ON REVERSE STRAND
            seqrev = record.reverse_complement()
            matchREV = re.search(str(seqrev.seq), mitogenome_seq, re.IGNORECASE)
            startREV = int(matchREV.start())
            endREV = int(len(ref.seq)) - int(matchREV.end())
            print(">" + record.id)
            print(("-" * startREV) + str(seqrev.seq) + ("-" * endREV))

        else:
            # GENE SEQUENCE IS ON DIRECT STRAND
            start = int(match.start())
            end = int(len(ref.seq)) - int(match.end())
            print(">" + record.id)
            print(("-" * start) + str(record_seq) + ("-" * end))


def main():
    args = get_args()
    mitogenome = args.mitofile
    genes = args.genesfile

    f = open(genes, 'r')
    # REMOVE EMPTY SPACES FROM FASTA HEADERS IN TEMPORARY FILES
    gclean = tempfile.NamedTemporaryFile(mode="w+", delete=False)
    txt = f.read().replace(' ', '')
    gclean.write(txt)
    gclean.close()
    f.close()

    if str(args.convertfile) == 'Y':
        # CONVERT AND SANITIZE HEADER FASTA NAMES IN TEMPORARY FILES
        gconverted = tempfile.NamedTemporaryFile(mode="w+",delete=False)
        for gene in SeqIO.parse(gclean.name, "fasta"):
            # SPLIT HEADER NAMES
            clean_name = str(gene.id).split(";")
            # SELECT 4TH FIELD
            new_name = clean_name[3]
            # SANITIZE NAME
            clean_new_name = re.sub("[\(\[].*?[\)\]]", "", new_name).replace("trn", "").replace("-", "_")
            gconverted.write(">" + clean_new_name + "\n" + str(gene.seq) + "\n")
        gconverted.close()
        # EXECUTE PROCESS
        process(mitogenome,gconverted.name)
        os.remove(gconverted.name)

    elif str(args.convertfile) == 'N':
        # EXEXUTE PROCESS
        process(mitogenome, gclean.name)

    else:
        print("ERROR: Type Y or N in --convertfile command")

    os.remove(gclean.name)



if __name__ == "__main__":
    main()
