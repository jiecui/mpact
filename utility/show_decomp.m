function show_decomp(s, P, fs, fig_name)
% SHOW_DECOMP Display the results of decomposition
%
% Syntax:
%
% Input(s):
%
% Output(s):
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2016 Richard J. Cui. Created: Sat 12/17/2016 10:37:41.537 AM
% $Revision: 0.1 $  $Date: Sat 12/17/2016 10:37:41.568 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

x       = hilbert(s);
N       = length(x);
nfreq   = 4 * N;
dbs     = 30; % dbs  range in dBs (optional, default is 25)
decf    = 1; % sub-sampling factor in time of the stft (must be integer)
w       = 05; % window width or window, must be a odd number

% Show fft spectrogram
% --------------------
[sx, t, f] = real_spec2(x, fs, nfreq, decf, w);   % spectrum without noise
show_tfd(sx, s, dbs, t, f);
set(gcf, 'Name', 'Spectrogram')

% show chirplet spectrogram
% ---------------------------------------
t_r = linspace(t(1), t(2), 2 * N); % time range
f_r = linspace(f(1), f(2), nfreq / 2 + 1); % frequency range
% esitmate WVD of chirplets directly
wrx = chirplet_spgrm(P, t_r, f_r, 'Method', 'direct');
show_tfd(wrx, s, dbs, t, f); % dB scale
set(gcf, 'Name', fig_name)

% use the explicit functions to calculate chirplet spectrogram
wig = chirplet_spgrm(P, t_r, f_r);
show_tfd(wig, s, dbs, t, f)
set(gcf, 'Name', fig_name)

end % function show_decomp

% [EOF]
