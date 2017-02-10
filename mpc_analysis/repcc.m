function [Pout,cc] = repcc(Pin,res)
% REPCC -- Rearrange chirplet parameters from decsend order of coherent
% coefficients.
% 
% syntax:
%   [Pout,cc] = repcc(Pin,res)
% 
% Inputs:
%   Pin     M by 5 matrix of chirplet parameters (see make_chirplets.m)
%   res     norm of the signal and the residuals for 1 to Q chirplets; 
%           could be used for a selection of stopping criterion
% 
% Outputs
%   Pout    rearranged Pin
%   cc      coherent coefficients whose size is one more than that of Pout.
% 
% The chirplet parameters (5 total) provided in M by 5 matrix. This
% function will rearrange them in a descend order according to the
% coherrent coefficients, which are defined as lambda = |A|/||res||, where
% A is the complex amplitude of the qth chirplet and ||res|| is the norm of
% the residue after extracting q chirplets from the signal x.
% 

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 28-Feb-2005
% $Revision: 0.1 $  $Date: 28-Feb-2005 11:04:57$

% Parse inputs and error check
error(nargchk(1,2,nargin));
if nargin ~= 2, error('Improper number of arguments for repcc!'); end

amp = abs(Pin(:,1));                % array of absolute amplitude
lambda = amp./res(2:end);
[cc,ix] = sort(lambda,'descend');   % coherent coefficients
Pout = Pin(ix,:);                        % chirplet parameters in new order

end