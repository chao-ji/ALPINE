function [result] = showresult(data, metric, scale)
%metric = 'correlation';
%scale = 'linear';
% metric function handel
mh = [];
if strcmp(metric, 'R-squared')	% r-squred
	mh = @rsquared;
elseif strcmp(metric, 'correlation')	% pearson correlation coefficient
	mh = @corr;
else
	error(fprintf(1, 'showresult: no such metric %s', metric));
end

% scale (linear or logarithm): transform function handle
th = [];
if strcmp(scale, 'log')
	th = @log10;
elseif strcmp(scale, 'linear')
	th = @(x) x;
else
	error(fprintf(1, 'showresult: no such transform function %s', scale));
end

% collecting results
Q = [];
for i = 1 : numel(data)
	% estimated protein quantity at Ith iteration
	Q = [Q; [data(i).Q, data(i).Qb, data(i).Qh, ...
	... % estiamted protein quantity of other metrics
	data(i).topn, data(i).ibaq, data(i).gm, data(i).meanint]];
end

fprintf(1, '%s under %s scale:\n', metric, scale);
% true protein quantity
tar = th(Q(:, 1));

result = cell(5 + numel(data(1).Qh) , 2);

for i = 0 : numel(data(1).Qh) 
%	fprintf(1, 'I = %d\t', i);
	result{i + 1, 1} = sprintf('I = %d', i);
end
result{numel(data(1).Qh) + 2, 1} = 'TopN';
result{numel(data(1).Qh) + 3, 1} = 'iBaq';
result{numel(data(1).Qh) + 4, 1} = 'Geomean';
result{numel(data(1).Qh) + 5, 1} = 'MeanInt';
%fprintf(1, 'TopN\tiBaq\tGeomean\tMeanInt\n');

for i = 2 : size(Q, 2)
	prd = th(Q(:, i));
%	fprintf(1, '%f\t', mh(prd, tar));
	result{i - 1, 2} = sprintf('%.4f', mh(prd, tar));
end

for i = 1 : size(result, 1)
	for j = 1 : size(result, 2)
		fprintf(1, '%s\t', result{i, j});
	end
	fprintf(1, '\n');
end
