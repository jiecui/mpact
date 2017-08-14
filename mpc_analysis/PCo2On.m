function P_on = PCo2On(N, P_co)
% PCO2ON convert chirplet parameters from Cohen to O'Neill format
%
% Syntax:
%   P_On = PCo2On(N, P_co)
% 
% Input(s):
%   N       - signal length
%   P_co    - chirplet parameters of Cohen format (see make_chirplets.m)
% 
% Output(s):
%   P_on    - chirplet parameters of O'Neill format
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
% $Revision: 0.2 $  $Date: Tue 05/30/2017 12:52:04.680 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

P_on = zeros(size(P_co));
% amplitude
P_on(:, 1) = P_co(:, 1);
% t center
P_on(:, 2) = P_co(:, 2)*N;
% f center
P_on(:, 3) = P_co(:, 3);
% chirp rate
P_on(:, 4) = P_co(:, 4)/N;
% duration
P_on(:, 5) = sqrt(N./P_co(:, 5)/2);

end % function PCo2On

% [EOF]
