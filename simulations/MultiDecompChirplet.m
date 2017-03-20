% SIM1_COMPOSITE_CHIRPLETS A simulation for different types of components
% 
% Syntax:
%   sim2_composite_chirplets
%
% Inputs:
%   None
%
% Outputs:
%   durkas      - The simulation signal without chirp
%   chirpsim    - The simulation signal with chirp = durkas+chirp
%
% Note:
%   This program produces the simulation signals used in the thesis to
%   compare with the signals used by Durka et al. References:
% 
% References:
%   [1] M. Akay and IEEE Engineering in Medicine and Biology Society, 
%       Time-frequency and wavelets in biomedical signal processing New
%       York: IEEE Press, 1998. pp. 305-406 
%   [2] P. J. Durka and K. J. Blinowska, "Analysis of EEG transients by
%       means of matching pursuit," Ann Biomed Eng, vol. 23, no. 5, pp.
%       608-611, Sept.1995.
%   [3]	J. Cui and W. Wong, "The adaptive chirplet transform and visual 
%       evoked potentials," IEEE Transactions on Biomedical Engineering,
%       vol. 53, pp. 1378-1384, Jul 2006.
% 
% See also .

% Copyright 2005-2016 Richard J. Cui. Created: Wed 03/20/2005  9:47:28.260 PM
% $Revision: 0.3 $  $Date: Fri 12/16/2016  3:11:09.454 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca


% =========================================================================
% create the signal
% =========================================================================
T = 512;        % the signal duration
fs = 1;         % sampling frequency

% ------
% wave I
% ------
waveI = zeros(1,T);

% component A
dt_a = 46; % even please, length of component A
tc_a = 90; % time-center of component A
fc_a = 1/dt_a; % frequency-center of A
s_a = -sin(2*pi*fc_a*(0:dt_a-1)/fs); % signal A

% component B, use the same as component A
tc_b = 422;
s_b1 = -(1/2 * sawtooth(2*pi*2*fc_a*(0:dt_a/2-1)/fs, 1/2) + 1/2);    
s_b2 = -s_b1;
s_b = [s_b1,s_b2]; % signal B

% component C
dt_c = 28; % length of C
fc_c = 0.4; % frequency center of C
tc_c = 256; % time center of C
s_c = gabor1d(T, fs, dt_c, fc_c, tc_c, 1, 0); % signal C

% composite wave I
% ----------------
waveI(round(tc_a-dt_a/2):round(tc_a-dt_a/2)+dt_a-1) = s_a;
waveI(round(tc_b-dt_a/2):round(tc_b-dt_a/2)+dt_a-1) = s_b;
waveI = waveI + s_c;

% -------
% wave II
% -------
% component D
dt_d = 2 * dt_c; % length of D
fc_d = 2 * fc_c / 3; % frequency center of D
tc_d = 256; % time center of D
s_d = gabor1d(T, fs, dt_d, fc_d, tc_d, 1, 0);
waveII = s_d;

% --------
% wave III
% --------
% component E - pulse
tc_e = 128;
s_e = gabor1d(T, fs, 0, 0, tc_e, 2, 0); % the pulse

% component F - sinusoidal wave
fc_f = .35;
A_f  = .2;
s_f  = A_f * sin(2*pi*fc_f * (0:T*fs-1)/fs);

% composite wave III
% ------------------
waveIII = s_e + s_f;
durkas = waveI + waveII + waveIII; % Durkas signal

% ----------------------------
% Secondly, add a chirp signal
% wave IV
% ----------------------------
% component G - chirplets
A_cp    = 6; % amplitude, total energy of the chirplet
tc_cp   = 350;  % time center
fc_cp   = .2 * 2 * pi / fs; % freqency center, unit rad
cr      = pi/T; % chirp rate
dt_cp   = 70;                % size of the chriplet
P   = [A_cp, tc_cp, fc_cp, cr, dt_cp];
cp  = make_chirplets(T, P);    % chirplets out is a complex column vector
waveIV = real(cp)';

% ------
% wave V
% ------
waveV = durkas + waveIV;
chirpsim = waveV; % composite signal

% =========================================================================
% draw the signals
% =========================================================================
% (1) draw Durka's signal and the composite signal
figure
subplot(211), plot(durkas), axis([0, T*fs-1, -2.5, 2.5]);
title('Durka''s signal without chirp');

subplot(212), plot(chirpsim), axis([0, T*fs-1, -2.5, 2.5]);
title('Composite signal with chirp component');

% (2) draw the structures of the signal
figure
subplot(414), plot(waveI), axis tight;
ylabel('waveI=A+B+C')

subplot(413), plot(waveII), axis tight;
ylabel('waveII=D')

subplot(412), plot(waveIII), axis tight;
ylabel('waveIII=E+F')

subplot(411), plot(waveIV), axis tight;
ylabel('waveIV=G')

% =========================================================================
% perform mp adaptive chirplet decomposition
% =========================================================================
% Common parameters
% -----------------
Q   = 7; % number of atoms desired
i0  = 1; % the first scale to roate the atoms
D   = 5; % decomposition depth = the higest scale
a   = 2; % the radix of scale
M   = T; % resolution for Newton-Raphson refinement
verbose = false; % show notes
mnits   = 5; % max number of iteration for refinement

% decompose waveV
% ---------------
P = mp_act_signal(waveV, Q, M, D, i0, a, 'ExpectMax', verbose, mnits);

% =========================================================================
% subroutines
% =========================================================================


% [EOF]
