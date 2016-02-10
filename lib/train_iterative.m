function [nets, data, mu, sigma] = train_iterative(data, par) 
% ---------------------------------------------------------------------
% Function: Trains an ensemble of neural networks for predicting peptide 
%	    respose rate iteratively. Each neural network is trained 
%	    separately on a bootstrapped sample.
%    Input: data -- struct array containing features and targets for 
%           peptides grouped by protein, data(i) is a struct that 
%           correspond to an individual protein which contains the 
%           feature (X), response rate (R), peak area (A) of peptides
%           and estimated quantity of protein (Q)
%   Output: nets -- A M by N cell array of neural networks, M = num of iterations
%	    N = num of bootstrapped sample
%           data -- the same as the input struct array, except that the
%           peptide response rate (R) and estimated quantity of protein (Q)
%           have been updated 
% ---------------------------------------------------------------------

nets = {};				% cell array of neural networks: rows -- number of iterations, cols -- number of bootstrapped samples
iterations = par.train.iterations;	% number of iterations
samples = par.train.samples;		% number of bootstrapped samples: Each sample contains the peptide ions from a SUBSET of PROTEINS
hidden_nodes = par.train.hidden_nodes;	% number of hidden nodes

% normalizing feature	
	X = [];
	for i = 1 : numel(data)
		X = [X; data(i).X];
	end
	[X, mu, sigma] = zscore(X);	size(X)
%	[X, proj] = pca_ret_var_noclass(X, 95);
	save('X', 'X');
	c = 1;
	for i = 1 : numel(data)
		data(i).X = X(c : (c + size(data(i).X, 1) - 1), :);
		c = c + size(data(i).X, 1);
	end

fprintf(1, '\n\nSTART ITERATIVE LEARNING OF RESPONSE RATE AND QUANTITY ESTIMATION.\n\n');

for i = 1 : iterations 

	fprintf(1, '\niteration = %d of %d\n', i, iterations);
	fprintf(1, '\nBootstrapping %d samples from %d training proteins\n', samples, numel(data));

	% During each iteration, use a subset of bootstrapped samples (the validation set) to estimate the optimal number of epoches
	[~, bootsam] = bootstrp(samples, [], data);
	fprintf(1, '\nUse validation set to determine the optimal number of epochs to train neural network\n');
	opt_epochs = getEpochs(data, bootsam, par);
	fprintf(1, '\nNum of epochs = %d\n', opt_epochs);

	X = {};	% cell array of features: each element corresponds to a boostrapped sample
	T = {}; % cell array of features: each element corresponds to a boostrapped sample
	fprintf(1, '\nStart training %d neural networks on %d bootstrapped samples\n', samples, samples);


	% For each j, a separate neural network will be using using the peptide ions from proteins in the jth bootstrapped sample
	for j = 1 : samples

%		bs = data(bootsam(:, i));
		bs = data(unique(bootsam(:, j)));
		X{j} = [];
		T{j} = [];
		for k = 1 : numel(bs)
			X{j} = [X{j}; bs(k).X];
			T{j} = [T{j}; log10(bs(k).R)];	%NOTE: the targets value (response rates) are log-transformed
		end		

		net = fitnet(hidden_nodes);
		net.trainFcn = 'trainbr';
		net.divideParam.trainRatio = 0.85;
		net.divideParam.val_ratio = 0.0;
		net.divideparam.testRatio = 0.15;
		net.trainParam.epochs = opt_epochs;
		fprintf(1, '\nStart training bootstrapped sample %d of %d\n', j, samples);
		fprintf(1, '\nNum of unique proteins in sample %d: %d\n', j, numel(bs));
		fprintf(1, '\nNum of peptide ions in sample %d: %d\n', j, size(X{j}, 1));
%		fprintf(1, '\nsize of X: %d by %d \n', size(X{i}, 1), size(X{i}, 2));
%		fprintf(1, '\nsize of T: %d by %d \n', numel(T{i}), 1);
		[net tr] = train(net, X{j}', T{j}');
%		ew = 1 ./ (10 .^ T{j}');
%		ew = ew / sum(ew);
%		[net tr] = train(net, X{j}', T{j}', [], [], ew);

		trainInd = tr.trainInd;
		testInd = tr.testInd;
		predTrain = net(X{j}(trainInd, :)');
		predTrain = reshape(predTrain, numel(predTrain), 1);
		predTest = net(X{j}(testInd, :)');
		predTest = reshape(predTest, numel(predTest), 1);
		fprintf(1, '\nIteration = %d, sample = %d: |training| = %d, trainPerf = %f and |test| = %d, testPerf = %f\n', i, j, numel(trainInd), corr(predTrain, T{j}(trainInd)), numel(testInd), corr(predTest, T{j}(testInd)));

		% 'net': the trained neural network for the ith iteration on the jth bootstrapped sample
		nets{i, j} = net;	
	end
	
	fprintf(1, '\nAn ensemble of %d neural networks have been trained\n', samples);
	
	% Do protein quantity estimation for each protein: j = 1 to length of 'data'
	for j = 1 : numel(data)
		X = data(j).X;
		predR = [];
		
		for k = 1 : samples
			net = nets{i, k};
			row = net(X');				%row: the predicted response rate (log-transformed) for peptide ions from the Jth protein using the Kth (k = 1 : samples) trained neural network in the Ith iteration
			row = reshape(row, 1, numel(row));	% 'row' converted to row vector
			predR = [predR; row];
		end

		if samples > 1
			predR = mean(predR, 1);			% average the predicted response rates
		end

		predR = 10 .^ predR; 				% transform the predicted log-transformed response rates by taking exponentiation 
		predR = reshape(predR, numel(predR), 1);	% 'predR' converted to col vector
		if any(predR == 0)
			error('train_iterative: Invalid value of predicted response rates');
		end

		Q_hat = geomean(data(j).A ./ predR);	
		data(j).Q(i) = Q_hat;				% compute the geometric mean of ratios (i.e. peak areas over response rates) as the updated protein quantity estimate
		data(j).R = data(j).A / data(j).Q(i);		% update the response rates
	end
end
