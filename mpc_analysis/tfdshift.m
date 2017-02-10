function out = tfdshift(in)
% tfdshift -- Shift the spectrum of a TFD by pi radians.
%
%  Usage
%    out = tfdshift(in)
%
%  Inputs
%    in   time-frequency distribution
%
%  Outputs
%    out  shifted time-frequency distribution

% Copyright 2005-2016 Richard J. Cui. Created: Tue 03/02/2005 2:46:52.278 PM
% $ Revision: 0.2 $  $ Date: Sun 12/11/2016 10:25:06.525 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

narginchk(1, 1);

N = size(in, 1);
M = ceil(N/2);
out = [in(M+1:N,:) ; in(1:M,:)];
