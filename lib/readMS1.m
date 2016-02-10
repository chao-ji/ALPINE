function [ms1] = readMS1(filename)
% ---------------------------------------------------------------------
% FUNCTION: read the file containing the peak area for each identified
%	    MS/MS spectra. Each line containes the Peak Area and the Scan
%           number (foreign key) corresponding to the MS/MS spectra
%   INPUT: filename -- file name
%  OUTPUT: a hash table with key (scan number) - value (a struct containing
%          scan number, peak area, begin and end scan number of the corresponding 
%          MS peak)
% ---------------------------------------------------------------------

if ~ischar(filename)
	error('readMS1:invalidUsage');
end

if ~exist(filename, 'file')
	error(sprintf('readMS1:FileNotFound %s\n', filename));
end

table = tdfread(filename);
ms1 = containers.Map;
for i = 1 : numel(table.FragScanNumber)
	element.scan = num2str(table.FragScanNumber(i));
	element.area = table.PeakArea(i);
	element.bnscan = table.PeakScanStart(i);
	element.edscan = table.PeakScanEnd(i);

	if ~ms1.isKey(element.scan)
		ms1(element.scan) = element;
	else
		error(sprintf('readMS1:DuplicateScan %s\n', element.scan));
	end
end
