% SIM2_X_CHIRPS A simulation of up- and down- chirps decomposition
%
% This simulation compares the results of EM and MLE refinement algorithm
% as the signal is embedded in noise.
% 
% Syntax:
%   sim_1_x_chirps
%
% Input(s):
%
% Output(s):
%
% Example:
% 
% Note:
%   For strong noise (SNR < 0 dB), MLE is more likely to fail than EM.
%
% See also .

% Copyright 2004-2017 Richard J. Cui. Created: Mon 09/27/2004  1:23:52.278 PM
% $ Revision: 1.5 $  $ Date: Mon 02/20/2017  9:17:09.916 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

%% ========================================================================
% Create and display the simulated signal with noise
% =========================================================================
d_snr = 0.0; % desired SNR
% synthesize the simulated signal
% -------------------------------
N   = 100; % signal size
t   = (0:N-1)'; % N seconds, sampling time = 1 sec
s1  = chirp(t, 0, N, .5); % up-chirp: 0 -> 0.5 Hz
s2  = flipud(s1); % down-chirp: 0.5 -> 0 Hz
s   = s1 + s2; % synthesized signal
[spn, ns] = add_noise(s, d_snr); % add noise to the signal (require communication toolbox)
signr = 10 * log10(var(s)/var(ns));
fprintf('Desired SNR = %.2f dB, estimated SNR = %.2f dB\n', d_snr, signr)

% display the signal
% ------------------
figure
% the synthesized signal
sh = subplot(414);
plot(spn), grid on, axis tight
axlm = axis(sh);
ay_max = max(floor(abs(axlm(3:4)) + .5));
ax_lmt = [axlm(1:2), [-1 1] * ay_max];
axis(sh, ax_lmt)
xlabel('Time (s)')
title(sprintf('S_1 + S_2 + noise (SNR = %.2f dB)', signr));
% signal s1
sh = subplot(411);
plot(s1), grid on, axis(sh, ax_lmt), title('Up-chirp S_1');
% signal s2
sh = subplot(412);
plot(s2), grid on, axis(sh, ax_lmt), title('Down-chirp S_2');
% s = s1 + s2
sh = subplot(413);
plot(s), grid on, axis(sh, ax_lmt), title('Synthesized signal');

%% ========================================================================
% perform mp adaptive chirplet decomposition
% =========================================================================
% Common parameters
% -----------------
Q   = 2; % number of atoms desired
i0  = 1; % the first scale to roate the atoms
D   = 5; % decomposition depth = the higest scale
a   = 2; % the radix of scale
M = 256; % resolution for Newton-Raphson refinement
verbose = false; % show notes
mnits   = 10; % max number of iteration for refinement
level   = 2; % difficulty of MLE (recommanded = 2)

% decompose spn
% -------------
P_em    = mp_act_signal(spn, Q, M, D, i0, a, 'ExpectMax', verbose, mnits);
P_mle   = mp_act_signal(spn, Q, M, D, i0, a, 'MaxLikeliEst', verbose, mnits, level);

% =========================================================================
% compare original and recontructed signals
% =========================================================================
% EM
fig_name = get_fig_name('ExpectMax');
comp_decomp(s, spn, P_em, fig_name)
% MLE
fig_name = get_fig_name('MaxLikeliEst');
comp_decomp(s, spn, P_mle, fig_name)

% =========================================================================
% subroutines
% =========================================================================

% [EOF]
