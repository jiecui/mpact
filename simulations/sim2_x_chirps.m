% SIM2_X_CHIRPS A simulation of up- and down- chirps decomposition
%
% This simulation compares the results of EM and MLE algorithm as the
% signal is embedded in strong noise
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
% $ Revision: 1.3 $  $ Date: Fri 02/10/2017  6:35:19.981 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

%% ========================================================================
% Create and display the simulated signal with noise
% =========================================================================
% synthesize the simulated signal
% -------------------------------
N   = 100; % signal size
t   = (0:N-1)'; % N seconds, sampling time = 1 sec
s1  = chirp(t, 0, N, .5); % up-chirp: 0 -> 0.5 Hz
s2  = flipud(s1); % down-chirp: 0.5 -> 0 Hz
s   = s1 + s2; % synthesized signal
d_snr = 0; % desired SNR
[spn, ns] = add_noise(s, d_snr); % add noise to the signal (require communication toolbox)
signr = 10 * log10(var(s)/var(ns));
fprintf('Desired SNR = %.2f dB, estimated SNR = %.2f dB\n', d_snr, signr)

% display the signal
% ------------------
figure
% the synthesized signal
subplot(311)
plot(spn), grid on, axis tight
ax_lmt = axis(gca);
title('Synthesized signal = S_1 + S_2 + noise');
% signal s1
subplot(312) 
plot(s1), grid on, axis(gca, ax_lmt), ylabel('Up-chirp S_1');
% signal s2
subplot(313) 
plot(s2), grid on, axis(gca, ax_lmt), ylabel('Down-chirp S_2');
xlabel('Time (s)')

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

% ========================================================================
% compare original and recontructed signals
% =========================================================================
comp_decomp(s, spn, P_em, 'ExpectMax')
comp_decomp(s, spn, P_mle, 'MaxLikeliEst')

% ========================================================================
% subroutines
% =========================================================================
function comp_decomp(s, spn, P, fig_name)
% compare original and recontructed signals

N = length(spn);
rs = real(make_chirplets(N, P)); % reconstructed signal

% compare original & reconstructed signals
% ----------------------------------------
figure('Name', fig_name)
subplot(311)
hs = plot(s);
axis tight, grid on;
hold on
hns = plot(spn);
hold off
legend([hs, hns], {'Clean' 'Noisy'})
ax_lim = axis(gca);
title 'Clean and noisy signal';

subplot(312)
hns = plot(spn);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(rs);
hold off
legend([hns, hrs], {'Noisy' 'Recon'})
title 'Noisy and reconstructed signal';

subplot(313), plot(spn-rs), axis(gca, ax_lim), grid on;
title 'Residuals = Noisy - Recon';

% compare clean & reconstructed signals
% ----------------------------------------
figure('Name', fig_name)
subplot(311)
hs = plot(s);
axis tight, grid on;
hold on
hns = plot(spn);
hold off
legend([hs, hns], {'Clean' 'Noisy'})
ax_lim = axis(gca);
title 'Clean and noisy signal';

subplot(312)
hs = plot(s);
axis(gca, ax_lim), grid on;
hold on
hrs = plot(rs);
hold off
legend([hs, hrs], {'Clean' 'Recon'})
title 'Clean and reconstructed signal';

subplot(313), plot(s - rs), axis(gca, ax_lim), grid on;
title 'Error = Clean - Recon';

end % function


% [EOF]
