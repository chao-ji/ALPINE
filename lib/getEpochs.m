function [opt_epochs] = getEpochs(data, bootsam, par)
% ---------------------------------------------------------------------
% Function: Determines optimal number of epochs to train neural network
%	    using a fraction of bootstrapped samples as validation set
%    Input: data -- struct array containing features and targets for 
%	    peptides grouped by protein, data(i) is a struct that 
%	    correspond to an individual protein which contains the 
%	    feature (X), response rate (R), peak area (A) of peptides
%	    and estimated quantity of protein (Q)
%	    bootsam -- M by N matrix, M = size of data, N = num of
%	    bootstrapped samples 
%	    par.train -- parameters for training
%   Output: Estimated optimal number of epcohs 
% ---------------------------------------------------------------------

val_ratio = par.train.val_ratio; % fraction of boostrapped samples to be including in the validation set
maxfail = par.train.maxfail;
hidden_nodes = par.train.hidden_nodes; % 6
epochs = par.train.maxepcohs; % 100 
samples = size(bootsam, 2); % number of bootstrapped samples

val_size = round(samples * val_ratio); % number of bootstrapped samples in the validation set
if val_size < 1
	opt_epochs = 100;
	fprintf(1, 'NOTE: getEpochs: Set number of epochs to default value 100.\n'); 
	return;
end

val_idx = sort(randsample(samples, val_size));

TR = [];
% each bootstrapped sample in the validation set is used to obtain an estimate of the number of epochs
% i.e. the number of epoches such that the performance achieved optimum in the test set
% 'opt_epochs' is calculated as the average of these estimates
for i = 1 : numel(val_idx)
	bs = data(unique(bootsam(:, val_idx(i))));
	X = [];
	T = [];
	for j = 1 : numel(bs)
		X = [X; bs(j).X];
		T = [T; log10(bs(j).R)];
	end
		
	net = fitnet(hidden_nodes); % neural network for non-linear regression
	net.trainFcn = 'trainbr';   % Bayesian regulation algorithm
	net.divideParam.trainRatio = 0.75;
	net.divideParam.valRatio = 0.0;
	net.divideparam.testRatio = 0.25;
	net.trainParam.epochs = epochs;
	[net tr] = train(net, X', T');
%	ew = 1 ./ (10 .^ T');
%	ew = ew / sum(ew);
%	[net tr] = train(net, X', T', [], [], 1 ./ (10 .^ T'));

	TR = [TR tr];
end

opt_epochs = [];
for i = 1 : val_size
	for j = 1 : numel(TR(i).tperf) - maxfail
		stop = [];	
		for k = j : j + maxfail - 1 
			stop = [stop, (TR(i).tperf(k) < TR(i).tperf(k + 1))];	
		end

		if all(stop) % stop if the squared error on the test set increases for 'maxfail' epochs in a row
			break;
		end
	end

	opt_epochs = [opt_epochs j];
end

if isempty(opt_epochs)
	error('getEpochs:Invalid estimate of optimal epochs\n');
end

opt_epochs = round(mean(opt_epochs));
