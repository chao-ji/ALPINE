function [data] = peppro2feature(protein, par)
% ---------------------------------------------------------------------
% FUNCTION: Converts a protein profile struct array into feature struct array
%    INPUT: protein -- a struct array containing protein profile information
%           (i.e. protein sequence, children peptides etc.) returned from function
%           buildProtProfile
%           par -- parameter struct
%   OUTPUT: a struct array containing features from each peptide (grouped by the
%           parent protein)
% ---------------------------------------------------------------------

neighbors = numel(par.feature.nr_range);
empty = repmat('-', 1, neighbors);
%%
%	Raw feature:
%	array of struct containing raw feature for EACH PROTEIN and its ID'ed peptides
RF = [];
%%
%	Feature:
%	vectorized values of raw feature
F = [];

for i = 1 : numel(protein) 

	rawfeature = struct('lNrAA', [], ... % neighboring amino acids on the left cleavage sites
			    'rNrAA', [], ... % neighboring amino acids on the right cleavage sites
					 ... % Example: ...NREY... the cleavage occurred on the 'R', 
					 ... % 'lNrAA' for the peptide 'EY...' will be 'NREY'
					 ... % 'rNrAA' for the peptide '...NR' will be 'NREY'
			     'lLoc', [], ... % amino acid indices of the left cleavage sites(K/R)
			     'rLoc', [], ... % amino acid indices of the right cleavage sites(K/R)
			   'pepseq', [], ... % the amino acid sequences of unique peptide ions to be used for quantifications 
			       'ch', [], ... % charge
			      'len', [], ... % protein sequence length
			     'area', []  ... % peak areas
			   ); 

	peptides = numel(protein(i).pep_for_quant);	% number of unique peptide ions to be used for quantification
	bn = [protein(i).pep_for_quant.bn]; % indices of N-term amino acid of unique peptide ions
	ed = [protein(i).pep_for_quant.ed]; % indices of C-term amino acid of unique peptide ions 
	% Example: protein = '...KPEPTIDER...', peptide = 'PEPTIDER'. 'bn' = index of the N-term 'P' and 'ed' = index of the C-term 'R'

	rawfeature.lLoc = bn - 1;		% indices of K/R
	rawfeature.rLoc = ed;			% indices of K/R
	rawfeature.len = numel(protein(i).sequence);
	rawfeature.lNrAA = cell(size(bn)); 
	rawfeature.rNrAA = cell(size(ed)); 
	rawfeature.pepseq = {protein(i).pep_for_quant.seq}; 
	rawfeature.ch = [protein(i).pep_for_quant.ch];

	for j = 1 : peptides 
		% extract left cleavage site neighboring residues
		neighbor = bn(j) - 1 + par.feature.nr_range;	% indices of the neighboring amino acids around the left cleavage site of the Jth peptide ion 
		valid = find(neighbor >= 1);	% make sure that the indices are valid (i.e. >= 1)
		nonempty = empty;
		nonempty(valid) = protein(i).sequence(neighbor(valid));
		rawfeature.lNrAA{j} = nonempty;

		% extract right cleavage site neitghboring residues
		neighbor = ed(j) + par.feature.nr_range;
		valid = find(neighbor <= rawfeature.len);
		nonempty = empty;
		nonempty(valid) = protein(i).sequence(neighbor(valid));
		rawfeature.rNrAA{j} = nonempty;
	end
	rawfeature.area = [protein(i).pep_for_quant.area];
	RF = [RF rawfeature];
end

for i = 1 : numel(RF)
i
	RF(i).feaPerAA = [predictBfactors(protein(i).sequence, 3, 5), VSL2B(protein(i).sequence), vihinen(protein(i).sequence)'];
end
save('RF', 'RF');
% vectorizing raw features (RF)

for i = 1 : numel(RF)
	peptides = numel(RF(i).pepseq);
	feature = struct(...
			 ...	% 1) Features describing the cleavage sites
			 'lNrAA', [], ...	% vectorized 'lNrAA' in raw features: containing indicator variables of the 
			 'rNrAA', [], ...	% vectorized 'rNrAA' in raw features
			  'lLoc', [], ...	
			  'rLoc', [], ...
			 ...	% 2) Features describing the amino acid composition of peptide ions
			     's', [], ...
			     'l', [], ...
			     'H', [], ...
			  'lens', [], ...
              'win0', [], ...
              'win5', [], ...
              'win10', [], ...
             'win15', [], ...
			 ...	% 3) Features describing the charge states
			    'ch', [], ...
			 ...	% Peak area			
			  'area', []);
	% left/right cleavage feature
	indicatorL = [];
	indicatorR = [];
	for j = 1 : peptides 
		indicatorL = [indicatorL; reshape(seq2indicator(RF(i).lNrAA{j})', 1, numel(par.feature.alphabet) * neighbors)];
		indicatorR = [indicatorR; reshape(seq2indicator(RF(i).rNrAA{j})', 1, numel(par.feature.alphabet) * neighbors)];
	end
	feature.lNrAA = indicatorL;	% indicator matrix of size numPeptides x 80 (i.e. 20 * neighbors)
	feature.rNrAA = indicatorR;	% indicator matrix of size numPeptides x 80 (i.e. 20 * neighbors)
	% normalize cleavage site indices by protein sequence length
	feature.lLoc = reshape(RF(i).lLoc, peptides, 1) / RF(i).len;	% matrix of size numPeptides x 1
	feature.rLoc = reshape(RF(i).rLoc, peptides, 1) / RF(i).len;	% matrix of size numPeptides x 1

	pepseq = RF(i).pepseq;
	% sequence feature
	% s: fixed-length vector representation of peptide sequences
	% l: peptide lengths
	% H: peptide entropies
	[feature.s, feature.l, feature.H] = seq2moment(pepseq, par.feature.alphabet, par.feature.moment_k);
	feature.ch = reshape(RF(i).ch, peptides, 1);		% matrix of size numPeptides x 1 
	feature.lens = repmat(RF(i).len, peptides, 1); 			% matrix of size numPeptides x 1

	p = RF(i).feaPerAA;
	for j = 1 : peptides
                win0 = [];
                b = RF(i).lLoc(j) + 1;
                e = RF(i).rLoc(j);
                b = b - 0;
                e = e + 0;
                if b < 1
                        b = 1;
                end
                if e > RF(i).len
                        e = RF(i).len;
                end
                win0 = [win0 mean(p(b:e, :))];
                feature.win0 = [feature.win0; win0];

                win5 = [];
                b = RF(i).lLoc(j) + 1;
                e = RF(i).rLoc(j);
                b = b - 5;
                e = e + 5;
                if b < 1
                        b = 1;
                end
                if e > RF(i).len
                        e = RF(i).len;
                end
                win5 = [win5 mean(p(b:e, :))];
                feature.win5 = [feature.win5; win5];

                win10 = [];
                b = RF(i).lLoc(j) + 1;
                e = RF(i).rLoc(j);
                b = b - 10;
                e = e + 10;
                if b < 1
                        b = 1;
                end
                if e > RF(i).len
                        e = RF(i).len;
                end
                win10 = [win10 mean(p(b:e, :))];
		feature.win10 = [feature.win10; win10];

                win15 = [];
                b = RF(i).lLoc(j) + 1;
                e = RF(i).rLoc(j);
                b = b - 15;
                e = e + 15;
                if b < 1
                        b = 1;
                end
                if e > RF(i).len
                        e = RF(i).len;
                end
                win15 = [win15 mean(p(b:e, :))];
                feature.win15 = [feature.win15; win15];
	end		

	feature.area = reshape(RF(i).area, peptides, 1);	% matrix of size numPeptides x 1
	F = [F feature];	
end

[data] = feature2data(F, RF, protein);
