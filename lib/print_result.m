function [] = print_result(data, outputfilename)
fid = fopen(outputfilename, 'w');

for i = 1 : numel(data)
	fprintf(fid, 'i = %d\n', i);
	fprintf(fid, 'Protein Accession: %s\t', data(i).accession);
	fprintf(fid, 'Estimated Quantity: %.2f\n', data(i).Q(end));
	fprintf(fid, '%d Peptide ions used for quantificatoin:\n', numel(data(i).pepions));
	for j = 1 : numel(data(i).pepions)
		fprintf(fid, '%s\t%f\n', data(i).pepions{j}, data(i).A(j));
	end
	fprintf(fid, '\n');
end
fclose(fid);
