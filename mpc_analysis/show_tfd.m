function show_tfd(tfd,x,dbs,t,f,fs,fnts,lin)
% SHOW_TFD -- show an image plot of a TFD with a linear or dB amplitude
% scale together with the original signal and its spectrum for a real signal.
%
%  Syntax:
%       show_tfddb(tfd,x,dbs,t,f,fs,fnts,lin)
%
%  Inputs:
%       tfd  time-frequency distribution
%       x    original signal, assume real signal
%       dbs  range in dBs (optional, default is 25)
%            scale
%       t    range of sampling times, unit seconds (optional)
%       f    range of frequency values, unit Hz (optional)
%       fs   sampling frequency (optional, default 1)
%       fnts font size of axis labels (optional)
%       lin  linear scale (optional, default lin =0, dB scale)

% Copyright 2005-2016 Richard J. Cui. Created: Tue 03/02/2005 2:46:52.278 PM
% $ Revision: 0.3 $  $ Date: Sun 12/11/2016 10:25:06.525 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% parse inputs
narginchk(1, 8);

if nargin < 8, lin = 0; end         % dB scale
if nargin < 7, fnts = 10; end
if nargin < 6, fs =1 ; end
if nargin < 5, f = [0.0 0.5]; end   % for real signal
if nargin < 4, t = [1 size(tfd,2)]; end
if nargin < 3, dbs = 25; end

if isempty(t), t = [1 size(tfd,2)]; end
if isempty(f), f = [0.0 0.5]; end
if isempty(fnts), fnts = 10; end

% find the portion of signal x
xa = floor(t(1)*fs)+1;
xb = floor(t(2)*fs)+1;
xseg = x(xa:xb);            % signal segment

% Find the fft of x

nfft = 1024;                % # points of fft
xfft = fft(x,nfft);         % fft of x 
if rem(nfft,2) == 0
    n = nfft/2+1;
else
    n = (nfft+1)/2;
end
famp = abs(xfft(1:n));
fa = floor(f(1)*nfft/fs)+1;
fb = floor(f(2)*nfft/fs)+1;
fampseg = famp(fa:fb);

% plot the image
if lin                              % linear amplitude scale
    mintfd = min(tfd(:));           % normailze between 0 and 1
    tfd = tfd - mintfd;
    maxtfd = max(tfd(:));
    tfd = tfd./maxtfd;
    clim = [0 1];
else                                % dB amplitude scale
    tfd = tfd./sum(sum(tfd));       % normalized the dynamic range to [0 1]
    tfd = 20*log10(abs(tfd)+eps);   % dBs
    tfd = tfd - max(max(tfd));
    clim = [-dbs,0];
end %if

disp_prest(tfd,t,f,xseg,fampseg,clim,fnts);

end