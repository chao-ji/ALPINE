function [vec] = seq2indicator(seq)
% ---------------------------------------------------------------------
% FUNCTION: Converts amino acid sequence to an indicator matrix
%    INPUT: seq -- a string from the 20 AA alphabet of length L 
%   OUTPUT: indicator matrix of size L x 20 
% ---------------------------------------------------------------------

vec = [];

for i = 1 : numel(seq)
	indicatorAA = zeros(1, 20);
	idx = aa2num(seq(i));
	if idx == 0	% '-'
		;
	elseif idx < 0	% invalid amino acid symbol
		error(sprintf('seq2indicator:Invalid amino acid symbol %s\n', seq(i)));
	else
		indicatorAA(idx) = 1;
	end
	vec = [vec; indicatorAA];
end

