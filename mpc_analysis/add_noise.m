function [y, n, snr_hat] = add_noise(x, dsnr)
% ADD_NOISE add guassian white noise to a signal
% 
%  Syntax:
%   y = add_noise(x, dsnr)
%
%  Inputs:
%   x       - signal
%   dsnr    - desired SNR in dB
% 
%  Outputs
%   y       - x + noise
%   n       - noise
%   snr_hat - estimated SNR in dB
% 
% Note:
%   The function requires Communication System Toolbox and Signal
%   Processing Toolbox.
% 
% See also .

% Copyright 2005-2017 Richard J. Cui. Created: Mon 03/02/2004  1:23:52.278 PM
% $ Revision: 0.5 $  $ Wed 02/22/2017 11:29:13.394 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

narginchk(1, 2);

% add noise
% ---------
y = awgn(x, dsnr, 'measured');
n = y - x;

% estimate SNR
% ------------
snr_hat = snr(x, n);

end % function

% [EOF]