function [y, n] = add_noise(x, dsnr)
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
% 
% Note:
%   The function requires Communication System Toolbox
% 
% See also .

% Copyright 2005-2016 Richard J. Cui. Created: Mon 03/02/2004  1:23:52.278 PM
% $ Revision: 0.4 $  $ Fri 12/16/2016 10:40:12.417 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

narginchk(1, 2);

y = awgn(x, dsnr, 'measured');
n = y - x;

end % function

% [EOF]