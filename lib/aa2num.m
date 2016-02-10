function [v] = aa2num(c)
% ---------------------------------------------------------------------
% FUNCTION: Converts amino acid symbols to index (between 1 and 20)
%    INPUT: Single captital letter representing an amino acid
%   OUTPUT: Corresponding index between 1 and 20 for a valid amino acid.
%           Or 0 if input is a '-'. See function 'peppro2feature'.
% --------------------------------------------------------------------- 
switch c
        case 'A'
                 v = 1;
        case 'R'
                 v = 2;
        case 'N'
                 v = 3;
        case 'D'
                 v = 4;
        case 'C'
                 v = 5;
        case 'E'
                 v = 6;
        case 'Q'
                 v = 7;
        case 'G'
                 v = 8;
        case 'H'
                 v = 9;
        case 'I'
                 v = 10;
        case 'L'
                 v = 11;
        case 'K'
                 v = 12;
        case 'M'
                 v = 13;
        case 'F'
                 v = 14;
        case 'P'
                 v = 15;
        case 'S'
                 v = 16;
        case 'T'
                 v = 17;
        case 'W'
                 v = 18;
        case 'Y'
                 v = 19;
        case 'V'
                 v = 20;
	case '-'
		 v = 0;
	otherwise
		 v = -1;
end

