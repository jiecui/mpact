function P_co = POn2Co(N, P_on)
% PON2CO convert chirplet parameters from O'Neill to Cohen format
%
% Syntax:
%   P_co = PCo2On(N, P_on)
% 
% Input(s):
%   N       - signal length
%   P_on    - chirplet parameters of O'Neill format (see make_chirplets.m)
% 
% Output(s):
%   P_co    - chirplet parameters of Cohen format
% 
% Example:
%
% Note:
%   In Cohen format, sampling frequency is normalized to N and the
%   frequency range is normalized to 2*\pi. In O'Neill format, sampling
%   frequency is normalized to 1 and frequency 2*\pi.
% 
% References:
%   [1] Cohen, L. (1995). Time-frequency analysis. Englewood Cliffs, N.J,
%       Prentice Hall PTR.
%   [2] O'Neill, J. C. and P. Flandrin (1998). Chirp hunting. Proceedings 
%       of the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis (Cat. No.98TH8380).
%
% See also make_chirplets.

% Copyright 2017 Richard J. Cui. Created: Mon 05/29/2017  4:31:00.890 PM
% $Revision: 0.1 $  $Date: Mon 05/29/2017  4:31:00.896 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

P_co = zeros(size(P_on));
% amplitude
P_co(:, 1) = P_on(:, 1);
% t center
P_co(:, 2) = P_on(:, 2)/N;
% f center
P_co(:, 3) = P_on(:, 3);
% chirp rate
P_co(:, 4) = P_on(:, 4)*N;
% duration
P_co(:, 5) = N./P_on(:, 5).^2/2;

end % function POn2Co

% [EOF]
