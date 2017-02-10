function show_tfddb(tfd,x,dbs,t,f,fnts)
% show_tfddb -- show an image plot of a TFD with a dB amplitude scale
% together with the original signal and its spectrum for a real signal.
%
%  Syntax:
%       show_tfddb(tfd,x,dbs,t,f,fnts)
%
%  Inputs:
%       tfd  time-frequency distribution
%       x    original signal, assume real signal
%       dbs  range in dBs (optional, default is 25)
%       t    vector of sampling times (optional)
%       f    vector of frequency values (optional)
%       fnts   font size of axis labels (optional)

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 26-Feb-2005
% $Revision: 0.2 $  $Date: 26-Feb-2005 22:32:54$

% parse inputs
error(nargchk(1, 6, nargin));

if nargin < 6, fnts = 10; end
if nargin < 5, f = [0.0 0.5]; end     % for real signal
if nargin < 4, t = [1 size(tfd,2)]; end
if nargin < 3, dbs = 25; end

if isempty(t)
    t = [1 size(tfd,2)];
end
if isempty(f)
    f = [0.0 0.5];
end

% Find the fft of x
N = length(x);
time = 1:N;

nfft = 1024;                % # points of fft
xfft = fft(x,nfft);         % fft of x 
if rem(nfft,2) == 0,
    n = nfft/2+1;
else
    n = (nfft+1)/2;
end
amp = abs(xfft(1:n));
f_range = linspace(0.0,0.5,n);  % small error here

% plot the image
tfd = tfd./sum(sum(tfd));       % normalized the dynamic range to [0 1]
tfd = 20*log10(abs(tfd)+eps);   % dBs
tfd = tfd - max(max(tfd));
disp_prest(tfd,t,f,x,amp,[-dbs,0],fnts);

end