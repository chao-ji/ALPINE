function [X, L, H] = seq2moment(seq, alphabet, K)
%%
%% Originally from Yong
%%
%% Input: 
%% 'seq': 1-by-num_pep cell array of peptide sequences 
%% 'alphabet': 1-by-num_AA array of amino acids
%% 'K': [0 1], moment of amino acid composition
%%
%% Intermediate variables:
%% 'L': 1-by-num_pep cell array of peptide lengths
%% 'indicatorM': 1-by-num_pep cell array of indicator matrices of peptide sequences
%%		 each indicator matrix is a pepLen-by-20 binary matrix
%% 'M0': 1-by-num_pep cell array of 1-by-20 matrices, stacking up the amino acid counts across the sequence
%% 'indicatorMNORM': same dimension as 'indicatorM', normalized
%%
%% 'M1': 1-by-num_pep cell array of 1-by-20 matrices

num_pep = length(seq);
num_AA = length(alphabet);
k = length(K);


%
L = cellfun(@length, seq, 'uniformoutput', 0);
indicatorM = cellfun(@(x) seq2indicator(x), seq, 'uniformoutput', 0);
% amino acid count 
M0 = cellfun(@(x) sum(x), indicatorM, 'uniformoutput', 0);
% normalized indicator matrix
indicatorMNORM = cellfun(@(x, y) x ./ repmat(y + (y == 0), size(x, 1), 1), ...
	indicatorM, M0, 'uniformoutput', 0);

if any(K > 0)
	M1 = cellfun(@(x) (1:size(x, 1)) * x, indicatorMNORM, 'uniformoutput', 0);
end


X = zeros(num_pep, num_AA * k);
for i = 1 : k
	if K(i) == 0
		X(:, (i - 1) * num_AA + (1 : num_AA)) = vertcat(M0{:});
	elseif K(i) == 1
		X(:, (i - 1) * num_AA + (1 : num_AA)) = vertcat(M1{:});
	else
		Mi = cellfun(@(x, y, z) sum((repmat((1 : z)', 1, num_AA) - repmat(y, z, 1))...
	 		.^ K(i) .* x) .^ (1 / K(i)), indicatorMNORM, M1, L,'uniformoutput', 0);
		X(:, (i - 1) * num_AA + (1 : num_AA)) = vertcat(Mi{:});
	end
end

L = vertcat(L{:});
H = myEntropy(vertcat(M0{:}), 2);
