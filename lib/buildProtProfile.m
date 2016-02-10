function [protein] = buildProtProfile(filename, pep_hits, protein_list)
% ---------------------------------------------------------------------
% FUNCTION: Builds profiles for proteins to be used for subsequent 
%           iterative learning of response rate and protein quantity
%           estimation
%    INPUT: filename -- FASTA file containing the amino acid sequences 
%	    of the proteins
%           pep_hits -- An array of struct containing peptide hits info
%           returned from the function 'parseResult'
%           protein_list -- A cell array of accession numbers of proteins 
%           to be included in the profile
%   OUTPUT: protein -- An array of struct containing the profiles of 
%           proteins (including the peptide hits to proteins)
% ---------------------------------------------------------------------
protein = [];
fasta = fastaread(filename);
prot_hash = containers.Map;

% start reading the protein sequences from FASTA file
for i = 1 : numel(fasta)
	if isempty(fasta(i).Header)
		error(sprintf('buildProtProfile:Protein %d Header is empty\n', i));
	end
	tmp = regexp(fasta(i).Header, ' ', 'split');
	accession = tmp{1};

	% include protein only if the protein accession appears in the list 'protein_list'
	if ~isempty(find(strcmp(accession, protein_list)))
		profile = struct('accession', accession,...		% protein accession
				  'sequence', fasta(i).Sequence,...	% amino acid sequence
				  'pep_hits', [],...			% peptide identifications mapped to this protein
				         'Q', [],...			% protein quantity to be estimated
			     'pep_for_quant', []);			% unique peptide ions to be used for quantification
		protein = [protein profile];
		prot_hash(accession) = numel(protein);
	end
end

if isempty(protein)
	error('buildProtProfile:The number of proteins to be used for iterative training is zero.\n');
end

% assign peptide hits to parent proteins
for i = 1 : numel(pep_hits)
	if numel(pep_hits(i).protein) ~= 1
		error(sprintf('buildProtProfile:Peptide hits (scan = %s) was not unique\n', pep_hits(i).spec));
	end

	accession = pep_hits(i).protein.id;
	
	if prot_hash.isKey(accession)
		idx = prot_hash(accession);
		protein(idx).pep_hits = [protein(idx).pep_hits pep_hits(i)];
	else
		fprintf(1, 'buildProteinProfile: peptide hit %s not found in protein %s \n', pep_hits(i).seq, accession);
	end
end

% map peptide sequences to protein sequences
for i = 1 : numel(protein)
				  
	pep_hits = protein(i).pep_hits;
	proseq = protein(i).sequence;

	% 'pepch_hash' hash: Key = 'peptide_charge', Value = indices in the 'pep_hits' corresponding
	% to the peptide ion 'peptide_charge'
	pepch_hash = containers.Map;
	for j = 1 : numel(pep_hits)
		str = [pep_hits(j).seq, '_', num2str(pep_hits(j).ch)];
		if ~pepch_hash.isKey(str)
			pepch_hash(str) = [j];
		else
			pepch_hash(str) = [pepch_hash(str) j];
		end
	end

	% 'pepch' contains the unique strings for peptide ions
	pepch = pepch_hash.keys;

	for j = 1 : numel(pepch)
		% 'peptide' struct containing 1) unique peptide ions (sequence, charge)
		%			      2) combined peak area
		%			      3) the begin and end scan number of the MS1 peak with the combined peak area
		%			      4) max and min score of the peptide hits corresponding to the unique peptide ion
	        peptide = struct('seq', [],...
	                          'ch', [],...
	                        'area', [],...
	                          'bn', [],...
	                          'ed', [],...
	                 'score_range', []);


		index = pepch_hash(pepch{j});
		peak_area = [pep_hits(index).area];

		% 'peak_range' hash: Key = 'bnscan_edscan', Value = peak_area
		% Note: each unique 'peptide_charge' may consists of multiple peptides hits, and each peptide hit is associated 
		% with an MS1 peak (defined by the begin_scan and end_scan in 'pep_hits'). The combined peak area is calculated 
		% as the sum of the peak areas of all the MS1 peaks. If multiple peptide hits correspond to the same MS1 peak
		% e.g. same bn_scan and ed_scan, then the maximum peak area of these peitde hits will be assigned to this MS1 peak
		% For example: 'PEPTIDE_3' corresponds to four peptide hits with (bn_scan, ed_scan, area):
		% (1000, 1200, 2000), (1000, 1200, 2000), (1250, 1300, 100) (1250, 1300, 50) 
		% The first two correspond to the same MS1 peak, and the peak area of the peak (1250, 1300) will be 100, so the
		% combined peak area equals 2000+100=2100

		peak_range = containers.Map;
		for k = 1 : numel(index)
			bnscan = num2str(pep_hits(index(k)).bnscan);
			edscan = num2str(pep_hits(index(k)).edscan);
			range = [bnscan, '_', edscan];

			if ~peak_range.isKey(range)
				peak_range(range) = [peak_area(k)];
			else
				peak_range(range) = [peak_range(range) peak_area(k)];
			end
		end

		allkey = peak_range.keys;
		for k = 1 : numel(allkey)
			key = allkey{k};
			peak_range(key) = max(peak_range(key));
		end
	
		tmp = regexp(pepch{j}, '_', 'split');
		peptide.seq = tmp{1};
		peptide.ch = str2num(tmp{2});
		peptide.area = sum(cell2mat(peak_range.values));

		pepseq = peptide.seq;
		bn = findstr(proseq, pepseq);

		% exclude peptides that map to multiple locations in the SAME PROTEIN 
		if numel(bn) ~= 1
			fprintf('buildProteinProfile: %s found in multiple locations in protein %s, skipped\n', pepseq, protein(i).accession);
			continue;
		end
		ed = bn + numel(pepseq) - 1;
		peptide.bn = bn;
		peptide.ed = ed;

		score = [pep_hits(index).ionscore];
		peptide.score_range = [max(score), min(score)];
		protein(i).pep_for_quant = [protein(i).pep_for_quant peptide];
	end	
end

%                                          %
% Keep proteins with non-zero peptide hits %
%protein = protein(find(cellfun(@(x) ~isempty(x.pep_hits), mat2cell(protein, [1], ones(1, numel(protein))))));
tmp = [];
for i = 1 : numel(protein)
	if ~isempty(protein(i).pep_for_quant)
		tmp = [tmp protein(i)];
	end
end
protein = tmp; clear tmp;
	
