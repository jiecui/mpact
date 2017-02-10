function [tfd,t,f] = spec2(x,fs,nfreq,decf,w,how)
% SPEC2 compute samples of the type II spectrogram.
%
%  Usage
%    [tfd,t,f] = spec2(x,fs,nfreq,decf,w,how)
%
%  Inputs
%    x      signal vector
%    fs     sampling frequency of x (optional, default is 1 sample/second)
%    nfreq  number of samples to compute in frequency (optional, default
%           is 256)
%    decf   sub-sampling factor in time of the stft (optional, default 
%           is 1, i.e. no sub-sampling)
%    w      if length(w)==1 then the window is a guassian with duration 'w'
%           (see chirplets.m), otherwise 'w' is the window.  'w' must have an
%           odd length to have a center point. (optional, default is a 
%           gaussian with a duration of 5)
%    how    if how='short' then the stft is computed for the same times
%           as the signal and 
%               size(tfd,2) = length(x)
%           if how='long' then the stft will be computed for times before
%           and after the times of the signal and
%               size(tfd,2) = length(x) + length(w) - 1
%           (optional, default is 'short')
%
%  Outputs
%    tfd  matrix containing the spectrogram of signal x (optional)
%    t    vector of sampling times (optional)
%    f    vector of frequency values (optional)
%
% If no output arguments are specified, then the spectrogram is
% displayed using ptfd(tfd, t, f).

% Copyright 2005-2016 Richard J. Cui. Created: Sat 03/02/2005 2:46:52.278 PM
% $ Revision: 0.2 $  $ Date: Sun 12/11/2016 10:25:06.525 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

narginchk(1, 6);
if (nargin == 1)
  [tfd, t, f] = stft2(x);
elseif (nargin == 2)
  [tfd, t, f] = stft2(x, fs);
elseif (nargin == 3)
  [tfd, t, f] = stft2(x, fs, nfreq);
elseif (nargin == 4)
  [tfd, t, f] = stft2(x, fs, nfreq, decf);
elseif (nargin == 5)
  [tfd, t, f] = stft2(x, fs, nfreq, decf, w);
elseif (nargin == 6)
  [tfd, t, f] = stft2(x, fs, nfreq, decf, w, how);
end

tfd = real(tfd.*conj(tfd));

if (nargout == 0)
  ptfd(tfd, t, f);
  clear tfd
end
