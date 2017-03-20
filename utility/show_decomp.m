function show_decomp(s, P, fs, fig_name)
% SHOW_DECOMP Display the results of decomposition on t-f plane
%
% Syntax:
%
% Input(s):
%   s           - the orignal (real) signal
%   P           - chirplet structure
%   fs          - sampling frequency
%   fig_name    - name of the figures
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

% Copyright 2016-2017 Richard J. Cui. Created: Sat 12/17/2016 10:37:41.537 AM
% $Revision: 0.4 $  $Date: Mon 03/20/2017 12:25:17.545 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

x       = hilbert(s);
N       = length(x);
nfreq   = 4 * N;
dbs     = 30; % dbs range in dBs (optional, default is 25)
decf    = 1; % sub-sampling factor in time of the stft (must be integer)
w       = 05; % window width or window, must be a odd number

% Show fft spectrogram
% --------------------
[sx, t, f] = real_spec2(x, fs, nfreq, decf, w);   % spectrum without noise
show_tfd(sx, s, dbs, t, f);
set(gcf, 'Name', 'Spectrogram')

% show chirplet spectrogram
% ---------------------------------------
rs = real(make_chirplets(N, P)); % reconstructed signal
t_r = linspace(t(1), t(2), 2 * N); % time range
f_r = linspace(f(1), f(2), nfreq / 2 + 1); % frequency range
% esitmate WVD of chirplets directly
wrx = chirplet_spgrm(P, t_r, f_r, 'Method', 'direct');
show_tfd(wrx, rs, dbs, t, f); % dB scale
set(gcf, 'Name', ['WVD - ', fig_name])

% use the explicit functions to calculate chirplet spectrogram
wig = chirplet_spgrm(P, t_r, f_r);
show_tfd(wig, rs, dbs, t, f)
set(gcf, 'Name', ['ACS - ', fig_name])

end % function show_decomp

% [EOF]
