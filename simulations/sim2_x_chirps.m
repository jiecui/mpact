% SIM2_X_CHIRPS A simulation of up- and down- chirps decomposition
%
% This simulation compares the results of MPEM and MLE algorithm to
% estimtate the signal embedded in noise.
% 
% Syntax:
%   sim2_x_chirps
%
% Input(s):
%
% Output(s):
%
% Example:
% 
% Note:
%   For strong noise, MLE is more likely to fail than MPEM.
%
% See also .

% Copyright 2004-2017 Richard J. Cui. Created: Mon 09/27/2004  1:23:52.278 PM
% $ Revision: 1.8 $  $ Date: Tue 03/07/2017 11:50:25.793 AM $
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
P1 = [10*exp(1i*0), N/2+1, pi/2,  pi/N, N/3]; % up-chirplet 0 -> pi
P2 = [10*exp(1i*0), N/2+1, pi/2, -pi/N, N/3]; % down-chirplet 0 -> -pi
s1 = real(make_chirplets(N, P1)); % the synthesized signal
s2 = real(make_chirplets(N, P2)); % the synthesized signal
s = s1+s2;
[spn, ns, signr] = add_noise(s, d_snr); % add noise to the signal (require communication toolbox)
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
title(sprintf('S_1 + S_2 + noise (SNR = %.2f dB)', d_snr));
% signal s1
sh = subplot(411);
plot(s1), grid on, axis(sh, ax_lmt), title('Up-chirp S_1');
% signal s2
sh = subplot(412);
plot(s2), grid on, axis(sh, ax_lmt), title('Down-chirp S_2');
% s = s1 + s2
sh = subplot(413);
plot(s), grid on, axis(sh, ax_lmt), title('clean = S_1 + S_2');

%% ========================================================================
% perform mp adaptive chirplet decomposition
% =========================================================================
% Common parameters
% -----------------
Q   = 2; % number of atoms desired
fs  = 1;

% decompose spn
% -------------
tests = hilbert(spn); % testing signal, complex
[~, P_mpem] = test_mpem_act(Q, tests);
[~, P_mle] = test_mle_act(Q, tests);

% =========================================================================
% compare original and recontructed signals
% =========================================================================
% t-f representation of clean signals
% -----------------------------------
P = [P1; P2];
show_decomp(s, P, fs, 'Clean signal')

% MPEM
% -----
p_mpem = table2array(P_mpem);
fig_name = get_fig_name('ExpectMax');
show_decomp(spn, p_mpem, fs, fig_name) % on t-f plane
comp_decomp(s, spn, p_mpem, fig_name) % on time domain

% MLE
% ---
p_mle = table2array(P_mle);
fig_name = get_fig_name('MaxLikeliEst');
show_decomp(spn, p_mle, fs, fig_name) % t-f plane
comp_decomp(s, spn, p_mle, fig_name) % on time domain

% compare square error
% --------------------
comp_sqerr(s, P_mpem, P_mle)

% =========================================================================
% subroutines
% =========================================================================

% [EOF]
