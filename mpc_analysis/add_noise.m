function [y, n, snr_hat] = add_noise(x, dsnr)
% ADD_NOISE add guassian white noise to a signal
% 
%  Syntax:
%   y = add_noise(x, dsnr)
%
%  Inputs:
%   x       - signal
%   dsnr    - desired SNR in dB. if 'Inf', no noise added.
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
% $ Revision: 0.6 $  $Date: Sun 02/26/2017  2:32:36.421 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

narginchk(1, 2);

% add noise
% ---------
if dsnr == Inf
    y = x;
else
    y = awgn(x, dsnr, 'measured');
end % if
n = y - x;

% estimate SNR
% ------------
snr_hat = snr(x, n);

end % function

% [EOF]