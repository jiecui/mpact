function [tfd, t, f] = real_spec2(x, fs, nfreq, decf, w, how)
% REAL_SPEC2 compute samples of the type II spectrogram for real signal
%
%  Syntax:
%    [tfd, t, f] = real_spec2(x, fs, nfreq, decf, w, how)
%
%  Inputs:
%    x      - signal vector
%    fs     - sampling frequency of x (optional, default is 1 sample/second)
%    nfreq  - number of samples to compute in frequency (optional, default
%             is 256)
%    decf   - sub-sampling factor in time of the stft (optional, default
%             is 1, i.e. no sub-sampling)
%    w      - if length(w)==1 then the window is a guassian with duration 'w'
%             (see chirplets.m), otherwise 'w' is the window.  'w' must have an
%             odd length to have a center point. (optional, default is a
%             gaussian with a duration of 5)
%    how    - if how = 'short' then the stft is computed for the same times
%             as the signal and 
%                   size(tfd,2) = length(x)
%             if how = 'long' then the stft will be computed for times before
%             and after the times of the signal and
%                   size(tfd,2) = length(x) + length(w) - 1
%             (optional, default is 'short')
%
%  Outputs
%    tfd    - matrix containing the spectrogram of signal x (optional)
%    t      - range of sampling times (optional)
%    f      - range of frequency values (optional)

% Copyright 2005-2016 Richard J. Cui. Created: Sat 02/26/2005 2:46:52.278 PM
% $ Revision: 0.5 $  $ Date: Sun 12/11/2016 10:25:06.525 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% Sparse inputs
narginchk(1, 6);

if nargin < 6, how = 'short'; end
if nargin < 5, w = 5; end
if nargin < 4, decf = 1; end
if nargin < 3, nfreq = 256; end
if nargin < 2, fs = 1; end

% get the spectrum
[tfd_spec2,t] = spec2(x,fs,nfreq,decf,w,how);   % could get f here

% get real part of the results
if rem(nfreq,2) == 0
    n = nfreq/2+1;
else
    n = (nfreq+1)/2;
end
tfd = tfd_spec2(nfreq-n+1:nfreq,:);
t = [t(1) t(end)];
f = [0 fs/2];

end % function

% [EOF]

