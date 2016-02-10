function [data] = eval_est(data, nets, mu, sigma)

for i = 1 : numel(data)

	data(i).Qh = zeros(1, size(nets, 1));
	data(i).Rh = zeros(numel(data(i).R), size(nets, 1));

	data(i).Qb = mean(data(i).A);
	data(i).Rb = data(i).A / data(i).Qb;

%	A = sort(data(i).A, 'descend');
%	if numel(A) >= 3
%		data(i).Qb = mean(A(1:3));
%	else
%		data(i).Qb = mean(A);
%	end
%	data(i).Rb = data(i).A / data(i).Qb;

	for j = 1 : size(nets, 1) % iterations
		predR = [];
		for k = 1 : size(nets, 2) % bootstrapped samples
			net = nets{j, k};

			% normalize test matrix using Mu and Sigma estiamted from the training matrix
			X = data(i).X;
			X = X - repmat(mu, size(X, 1), 1);
			X = X ./ repmat(sigma, size(X, 1), 1);
	
			row = net(X');
			row = reshape(row, 1, numel(row));
			predR = [predR; row];
		end

		if size(nets, 2) > 1
			predR = mean(predR, 1);
		end

		predR = 10 .^ predR;
		predR = reshape(predR, numel(predR), 1);
		if any(predR == 0)
			error('eval_est: Invalid value of predicted response rates');
		end
		Qh = geomean(data(i).A ./ predR);
		data(i).Qh(j) = Qh;
		data(i).Rh(:, j) = data(i).A / data(i).Qh(j);
	end
end
