function [P1,cc1] = cohstr(P0,cc0,tlb,tub)
% COHSTR -- Find coherent structure whose coherent coefficients are between
% selected thresholds.
%
% syntax:
%     [P1,cc1] = cohstr(P0,cc0,tlb,tub)
%
% Inputs:
%     P0      input M by 5 matrix of chirplet parameters (see make_chirplets.m)
%     cc0     input coherent coefficients
%     tlb     lower bound of the threshold
%     tub     upper bound of the threshold
%
% Outputs:
%     P1      output chirplet parameters
%     cc1     output coherent coefficients
%

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 28-Feb-2005
% $Revision: 0.2 $  $Date: 01-Mar-2005 10:22:42$

% Parse inputs and error check
error(nargchk(1,4,nargin));
if nargin ~= 4, error('Improper number of arguments for cohstr!'); end

idx = cc0 >= tlb & cc0 <= tub;
Q = sum(idx);
if Q == 0
    fprintf('No cc is between %g and %g\n',tlb,tub);
    P1 = [];
    cc1 = [];
else
    P1 = P0(idx,:);
    cc1 = cc0(idx);
end %if

end