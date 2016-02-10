function [r2] = rsquared(prd, tar)
if numel(prd) ~= numel(tar)
	error(fprintf(1, 'rsquared: sizes of predicted and target values do not agree.\n'));
end     
len = numel(prd);

% single-variable linear regression
% b: slope 
b = sum((prd - mean(prd)) .* (tar - mean(tar))) / sum((prd - mean(prd)) .^ 2);
% a: intercept
a = mean(tar) - mean(prd) * b;

% compute predicted value through linear regression
pred = prd * b + a;

% sse: sum of squared error
% sst: sum of squared total
sse = 0;
sst = 0;
for i = 1 : len
	sse = sse + (pred(i) - tar(i)) ^ 2;
	sst = sst + (tar(i) - mean(tar)) ^ 2;
end
r2 = 1 - sse / sst;

