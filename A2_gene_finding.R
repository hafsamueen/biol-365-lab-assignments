# BIOL 365 Assign2 Gene Annotation using DECIPHER
# < Hafsa Mueen >, 20947674 
# This assignment uses DECIPHER and dependent R packages 
R.version.string
date()

# Task 1 - Loading and setting up the library and work directory used for assignment 2
library(DECIPHER)            # loading in DECIPHER 
setwd("~/Documents/BIOL365") # setting up work directory 

# Loading and reading the A2contig fasta file for annotation
A2contig = readDNAStringSet("~/Documents/BIOL365/A2contig.fasta", format="fasta") # reading FASTA file
A2contig # checking to see if the data loaded properly
names(A2contig) # reading the FASTA file header for the species name, Leptospira alexanderi

# Task 2 - reads FASTA file to predict the start and stop positions of protien coding genes in the given genome
orfs = FindGenes(A2contig, showPlot = TRUE, allScores = TRUE) # predicts start and stop positions

# subsetting genome to only include ORF's that are predicted as genes
genes = orfs[orfs[, "Gene"]==1,]    # genes are flagged as 1
dna = ExtractGenes(genes, A2contig) # extracting genes from the genome

# looking at the distribution of predicted start codons
table(subseq(dna[-1], 1, 3)) # distribution of predicted start codons
w = which(!subseq(dna, 1, 3) %in% c("ATG", "GTG", "TTG")) # checking which sequences have non-canonical bacterial start codons 
dna = dna[-w]     # removing said sequences from DNA dataset
genes = genes[-w] # removing said sequences from genes dataset


# converting sets to data frames + writing csv files to check full seq's.
dna_df = data.frame(width = width(dna), seq = as.character(dna))   
write.csv(dna_df, file='~/Documents/BIOL365/dna_seq.csv') # to find shortest and longest gene
print(genes) # to find start site locations and orientations

# translating genes to look at predicted protien sequences
amino_acids = ExtractGenes(genes, A2contig, type="AAStringSet") # extracting as amino acids instead of DNA
amino_acids
aa_df = data.frame(width = width(amino_acids), seq = as.character(amino_acids)) # creating dataframe
write.csv(aa_df, file='~/Documents/BIOL365/aa_seq.csv') # to find corresponding protien sequences of genes

# Task 3 - searching the contig for any non-coding genes
data("NonCodingRNA_Bacteria")
x = NonCodingRNA_Bacteria # storing dataset
names(x) # checking available RNA types

rnas = FindNonCoding(x, A2contig) # finding non-coding RNA in genome
rnas                              # printing predicted non-coding RNA type

# extracting and matching identified RNA type with annotations
annotations = attr(rnas, "annotations")
m = match(rnas[, "Gene"], annotations)
sort(table(names(annotations)[m])) # printing the different RNA types found

rna_genes = ExtractGenes(rnas, A2contig, type="RNAStringSet")
rna_genes # printing extracted RNA seq.

# integrating non coding RNA into protien coding gene predictions - doesnt work..
combined_genes = FindGenes(A2contig, includeGenes = rnas, allScores = TRUE)
combined_genes # printing coding and non coding genes together
