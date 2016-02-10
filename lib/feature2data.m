function [data] = feature2data(F, RF, protein)
% ---------------------------------------------------------------------
% FUNCTION: Combines the features calculated for each peptide (amino acids
%           around the cleavage sites, amino acid composition of the 
%           peptide sequence, and the charge etc.)
%	    Caculates the initial quantity estimates for each protein (e.g.
%           average of peak areas of peptide ions per protein), and the 
%           initial estimates of peptide response rates
%    INPUT: F -- An array of struct containing the features for each peptide (grouped by protein)
%           RF -- An array of struct containing the raw features for each peptide (grouped by protein)
%	    Both 'F' and 'RF' are returned from the function 'peppro2feature'
%	    protein -- An array of struct containing the protein profile information
%   OUTPUT: data -- An array of struct containing the feature vectors.
%           Each data(i) corresponds to a protein: data(i).X: the matrix
%           containing features (num of rows = num of peptides, num of cols 
%           = num of features), data(i).A: peak area, data(i).Q:
%           initial quantity estimation, and data(i).R: the initial 
%           response rate estimation.
% ---------------------------------------------------------------------

%save('tmp', 'F', 'RF', 'protein');
for i = 1 : length(F) 
	data(i).X = [F(i).lNrAA F(i).lLoc F(i).s F(i).l F(i).H F(i).ch F(i).win0 F(i).win5 F(i).win10 F(i).win15 F(i).rNrAA F(i).rLoc F(i).lens];
	% 'mean(F(i).A)' : average peak area as initial estimation of protein quantity
	data(i).A = F(i).area;
	data(i).Q = mean(F(i).area);
	data(i).R = F(i).area / mean(F(i).area);
	data(i).accession = protein(i).accession;
	data(i).sequence = protein(i).sequence;
	data(i).pepions = cellfun(@(x) [x.seq, '_', num2str(x.ch)], num2cell(protein(i).pep_for_quant), 'UniformOutput', 0); 
end
