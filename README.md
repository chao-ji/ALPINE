# ALPINE
An EM-like algorithm that estimates protein abundance and peptide response rate simultaneously

Absolute and Label-free Protein quantification based on ion INtEnsity and predicted peptide response rate
School of Informatics and Computing, Indiana University

Contact: Chao Ji, jic@indiana.edu 

1. Input: 
	a. file containing peptide peak area and associated MS2 scan
number (A10-12023_SICstats.txt)
	b. file containing peptide identification results (F002994.xml)
	c. file containing accession numbers of proteins to be included in the
iterative algorithm (protein_list.txt)
	d. protein sequence file (pompep.fasta)

2. Usage:
	In the matlab environment execute pipeline.m
	>> pipeline

3. Ouput:
	The output file contains estimated protein abundances and the peptides
(with peak area) used for quantification

