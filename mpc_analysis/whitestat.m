function wm = whitestat(x)
% WHITESTAT - A staticstic of whiteness measure
%
% Syntax:
%   wm = whitestat(x)
%
% Inputs:
%   x           stochastic sequence, assume i.i.d. If it is complex, just
%               consider the real part.
%   
% Outputs:
%   wm          whiteness measure
%
% References:
%   K. Drouiche, "A new test for whiteness", IEEE Trans. Sig. Proc., Vol.
%   48, No. 7, July 2000.

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 13-Apr-2005
% $Revision: 0.2$  $Date: 13-Apr-2005 21:14:36$

% parameter settings
x = real(x);        % just get real part
euler = 0.57721;    % Euler's constant
N = length(x);

% find the covariance
x = x(:);       % force column vector
r0 = cov(x);    % autocovariance

% find periodogram
IN = abs(fft(x)).^2/2/pi/N;
% equivalently, we can use Matlab's function PERIODOGRAM
% IN = periodogram(ns,[],N,1,'twosided');

% outputs
wm = log(r0/2/pi)-sum(log(IN))/N-euler;
% equivalently
% wm = log(r0)-sum(log(IN))/N-euler;

end