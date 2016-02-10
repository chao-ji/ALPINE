function [pep_hits] = parseResult(filename, par)
% ---------------------------------------------------------------------
% FUNCTION: Parses the peptide identification result file in XML format
%    INPUT: filename -- The name of the peptide identification
%           result file
%	    par -- Struct containing parameters for filtering input data,
%           extracting features, and training the neural network
%   OUTPUT: An array of struct containing peptide hits info extracted
%           from input file 'filename'
%     NOTE: Only UNIQUE peptides (i.e. mapped to one parent protein) will 
%           be included
% ---------------------------------------------------------------------
if ~ischar(filename)
	error('parseResult:invalidUsage');
end

if ~exist(filename, 'file')
	error(sprintf('parseResult:FileNotFound %s\n', filename));
end

strct = parseXML(filename);
ms1info = readMS1(par.input.ms1file);
pep_hits = [];

for i = 1 : numel(strct.Children(2).Children)
i
	if strcmp(strct.Children(2).Children(i).Name, 'spectrum_query') && ~isempty(strct.Children(2).Children(i).Children)
		hit = struct('seq', [],...
			    'spec', [],...
                            'scan', [],...
                              'ch', [],...
                        'ionscore', [],...
                      'misscleave', [],...
                           'modif', [],...
                         'protein', [],...
                            'area', [],...
                          'bnscan', [],...
                          'edscan', []);
		
		for j = 1 : numel(strct.Children(2).Children(i).Attributes)
			name = strct.Children(2).Children(i).Attributes(j).Name;
			value = strct.Children(2).Children(i).Attributes(j).Value;
			if strcmp(name, 'spectrum')
				tmp = regexp(value, ' ', 'split');
				hit.spec = tmp{1};
				tmp = regexp(value, '\.', 'split');
				hit.scan = tmp{2};
			elseif strcmp(name, 'assumed_charge')
				hit.ch = str2num(value);
			end
		end

		for j = 1 : numel(strct.Children(2).Children(i).Children(2).Children(2).Attributes)
			name = strct.Children(2).Children(i).Children(2).Children(2).Attributes(j).Name;
			value = strct.Children(2).Children(i).Children(2).Children(2).Attributes(j).Value;

 			if strcmp(name, 'hit_rank')
				if ~strcmp(value, '1')	% select only top hits
					error(sprintf('parseXML:NoTopHits %s\n', hit.spec));
				end
			elseif strcmp(name, 'num_missed_cleavages')
				hit.misscleave = str2num(value);
			elseif strcmp(name, 'peptide')
				hit.seq = value;
			elseif strcmp(name, 'protein')
				hit.protein = [hit.protein struct('id', value)];
			end
		end

		for j = 1 : numel(strct.Children(2).Children(i).Children(2).Children(2).Children)
			tag = strct.Children(2).Children(i).Children(2).Children(2).Children(j);
			if strcmp(tag.Name, 'search_score')
				if ~isempty(find(strcmp({tag.Attributes.Value}, 'ionscore')))
					idx = find(strcmp({tag.Attributes.Name}, 'value'));
					hit.ionscore = str2num(tag.Attributes(idx).Value);
				end
			elseif strcmp(tag.Name, 'alternative_protein')
				idx = find(strcmp({tag.Attributes.Name}, 'protein'));
				hit.protein = [hit.protein struct('id', tag.Attributes(idx).Value)];
			elseif strcmp(tag.Name, 'modification_info')
				for k = 1 : numel(tag.Children)
					if strcmp(tag.Children(k).Name, 'mod_aminoacid_mass')
						mod = struct('modpos', [], 'modmass', []);

						idx = find(strcmp({tag.Children(k).Attributes.Name}, 'position'));
						mod.modpos = str2num(tag.Children(k).Attributes(idx).Value);
						idx = find(strcmp({tag.Children(k).Attributes.Name}, 'mass')); 
						mod.modmass = str2num(tag.Children(k).Attributes(idx).Value);

						hit.modif = [hit.modif mod];
					end
				end	

			end
		end
		
		if ms1info.isKey(hit.scan)
			ms1 = ms1info(hit.scan);
			hit.area = ms1.area;
			hit.bnscan = ms1.bnscan;
			hit.edscan = ms1.edscan;
		else
			error(sprintf('parseResult:NoMS1Info %s\n', hit.scan));
		end

		% include if 1) the number of parent protein = 1
		%            2) ion score >= cutoff
		if hit.ionscore >= par.input.cutoff && numel(hit.protein) == 1
			% check the parameter setting w.r.t if modified peptides will be included
			if par.input.no_modif
				% include if the peptide does not contain modification
				if isempty(hit.modif)
					pep_hits = [pep_hits hit];
				end
			else
				% include no matter if the peptide contains modification or not
				pep_hits = [pep_hits hit];
			end
		end
	end
end

