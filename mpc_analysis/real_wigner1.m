function [wig,t,f] = real_wigner1(x,fs,nfreq)
% REAL_WIGNER1 -- Compute samples of the (type I) Wigner distribution
% for real signal.
%
%  Syntax:
%    [wig,t,f] = real_wigner1(x,fs,nfreq)
%
%  Inputs
%    x      signal vector.  wigner1 assumes that x is sampled at the Nyquist
%           rate and uses sinc interpolation to oversample by a factor of 2.
%           If there are two columns then a cross Wigner distribution is
%           computed.
%    fs     sampling frequency of x (optional, default is 1 sample/second)
%    nfreq  number of samples to compute in the frequency direction, must
%           be at least 2*length(x) (optional, defaults to 2*length(x))
%
%  Outputs
%    wig    matrix containing the Wigner distribution of signal x.  If x has
%           length N, then tfd will be 2N by 2N. (optional)
%    t      vector of sampling times (optional)
%    f      vector of frequency values (optional)

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 27-Feb-2005
% $Revision: 0.1 $  $Date: 27-Feb-2005 10:39:36$

% parse inputs
narginchk(1,3);

N = length(x);      % signal size
if nargin < 3, nfreq = 2*N; end
if nargin < 2, fs = 1; end
if nfreq < 2*N,error('nfreq must be at least 2*length(x)'); end

%% Main body
% create the temporal acf for positive tau
acf = lacf1(x);     % type I local acf
NN = 2*N;           % width of acf

% negative tau
acf = [acf ; zeros(nfreq+1-NN, NN) ; conj(flipud(acf(2:NN/2,:)))];
% WVD
wig_wigner1 = real(fft(acf))/NN;
wig_wigner1 = tfdshift(wig_wigner1);

t = 1/(2*fs) * (0:NN-1);

% get real part of the results
if rem(nfreq,2) == 0
    n = nfreq/2+1;
else
    n = (nfreq+1)/2;
end
wig = wig_wigner1(nfreq-n+1:nfreq,:);
t = [t(1) t(end)];
f = [0 fs/2];

end % function

% [EOF]