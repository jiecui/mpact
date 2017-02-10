function g = gabor1d(T, fs, dt, fc, tc, A, phi)
% GABOR1D constructs 1-D Gabor wave
% 
% Syntax:
% 
% Input(s):
%   T       - signal duration
%   fs      - sampling frequency
%   dt      - time-duration/spread of the Gabor
%   fc      - frequency center of the Gabor
%   tc      - time center of the Gabor
%   A       - amplitude of the Gabor (real number)
%   phi     - the phase of the Gabor wave (real number, in rad)
%
% Output(s):
%   g       - the Gabor
% 
% Note:
%   All values in physical units: seconds, Hz etc.
%   (1) Sinusoid: dt = T
%   (2) Pulse: dt = 0
%   (3) Gaussian: fc = 0
% 
% Example:
%   g = gabor1d(512, 1, 28, .4, 256, 1, 0)
% 
% See also .

% Copyright 2016 Richard J. Cui. Created: Fri 12/16/2016  4:35:53.730 PM
% $ Revision: 0.1 $  $ Date: Fri 12/16/2016  4:35:53.730 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% =========================================================================
% normalize parameters to samples (points):
% =========================================================================
t = 0:(T * fs - 1);
p = tc * fs;
w = dt * fs;
f = fc / fs;

% =========================================================================
% construct the signa;
% =========================================================================
if dt == T % sinusoid
    g = A * cos(2 * pi * f .* t + phi);
elseif dt == 0 % pulse
    g = zeros(size(t));
    g(round(p)) = A;
elseif fc == 0 % a Gaussian
    g= A * exp(-((t - p) / w) .^ 2);
else
    g = A * exp(-((t - p) / w) .^ 2) .* cos(2 * pi * f .* t + phi);
end % if

end % function gabor1d

% [EOF]