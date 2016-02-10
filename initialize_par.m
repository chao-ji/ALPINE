% parameters for processing input file 
par.input.ms1file = 'A10-12023_SICstats.txt';	% input file containing peptide peak area and associated MS2 scan number
par.input.pepIDfile = 'F002994.xml';	% Mascot file containing peptide identification results
par.input.fastafile = 'pompep.fasta';	% protein sequence file
par.input.protein_list = 'protein_list.txt';	% protein accession numbers to be included in the iterative algorithm
par.input.cutoff = 0;	% MASCOT ionscore cutoff
par.input.no_modif = 1;	% 1: not including peptide with modification, 0: including peptide with modification

% parameters for extracting features from peptide and protein sequence (DO NOT change)
par.feature.nr_range = -1:2;
par.feature.alphabet = 'ARNDCEQGHILKMFPSTWYV';
par.feature.moment_k = [0 1];

% parameters for the iterative algorithm
par.train.iterations = 3;	% number of the iterations 
par.train.samples = 30;	% number of bootstraping samples
par.train.hidden_nodes = 10; 	% number of hidden neurons
par.train.val_ratio = 0.3;	% fraction of boostraping samples to be including in the validation set
par.train.maxfail = 6;	
par.train.maxepcohs = 150;
par.outputfilename = 'output.txt';	% output file containing estimated protein abundance and peptides used to quantify proteins
