function H = myEntropy(p,dim)

% entropy computation
% dim is the dimension for probability normalization
pSum = sum(p,dim);
d = ones(1,numel(size(p)));
d(dim) = size(p,dim);
p = p./repmat(pSum,d);
p = p.*log2(1./p);
p(~isfinite(p)) = 0;
H = sum(p,dim);
