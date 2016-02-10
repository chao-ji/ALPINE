addpath('~/newALPINE/lib');
addpath(genpath('~/matlab'));

initialize_par;

[pep_hits] = parseResult(par.input.pepIDfile, par);
save('pep_hits', 'pep_hits');
