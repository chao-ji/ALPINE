function [data] = normdata(data)

X = [];

for i = 1 : numel(data)
	X = [X; data(i).X];
end
X = zscore(X);

count = 1;
for i = 1 : numel(data)
	data(i).X = X([count : count + size(data(i).X, 1) - 1], :);
	count = count + size(data(i).X, 1);
end
