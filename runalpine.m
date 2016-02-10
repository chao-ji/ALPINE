addpath('~/ALPINE/lib');
initialize_par;

load pep_hits_filtered.mat;
list = importdata(par.input.protein_list);
protein = buildProtProfile(par.input.fastafile, pep_hits, unique(list));
save('protein', 'protein');
[data] = peppro2feature(protein, par);
save('data', 'data');

[nets, data, mu, sigma] = train_iterative(data, par);
save('nets', 'nets', 'data', 'mu', 'sigma');

load trueabd.mat;
list = allquan.keys;
refpro = buildProtProfile(par.input.fastafile, pep_hits, list);
data_ref = peppro2feature(refpro, par);
for i = 1 : numel(data_ref)
	data_ref(i).Q = allquan(data_ref(i).accession);
	data_ref(i).R = data_ref(i).A / data_ref(i).Q;
end
data_ref = eval_est(data_ref, nets, mu, sigma);
data_ref = addothermetric(data_ref);
showresult(data_ref, 'correlation', 'log');
save('data_ref', 'data_ref');


