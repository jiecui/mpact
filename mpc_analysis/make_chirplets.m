function x = make_chirplets(N, P, dflag)
% MAKE_CHIRPLETS construct a signal that is a sum of chirplets
%
%  sytax:
%   x = chirplets(N, P)
%
%  inputs:
%   N       - length of signal 
%   P       - matrix of parameters [A t f cr d; ...],
%               A:  complex amplitude = |A|e^{j\phi}
%               t:  time center (sample or unit sampling time)
%               f:  frequency center (rad)
%               cr: chirprate (rad/sample)
%               d:  chirplet duration (sample)
%             (optional, default is [1 N/2+1 0 0 sqrt(N/4/pi)])
%
%   dflag   - consider the discretization effect (Optional: default 0)
% 
%  outputs:
%   x       - constructed signal
%
% Note that d is the standard deviation of the guassian, d=sqrt(N/4/pi)
% gives an atom with a circular Wigner distribution, and 2*sqrt(2)*d is the
% Rayleigh limit.
%
% Assume real signals & sampling frequency = 1 Hz. Note that Nyquist = 1
% for Matlab conventions. Use rad for frequency unit.
%
%  Examples
%    N = 128; x = make_chirplets(N, [1 N/2+1 0 0 sqrt(N/4/pi)]);
%    x = make_chirplets(128, [1 55 0 2*pi/128 12 ; 1 75 0 2*pi/128 12]);

% Copyright 2005-2017 Richard J. Cui. Created: Thu 03/02/2005 10:33:22.814 AM
% $ Revision: 0.4 $  $ Date: Tue 02/21/2017  1:53:17.867 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com


narginchk(1,3);

if ~exist('dflag', 'var')
    dflag = 0;      % consideration of discretization effect (problem?)
end % if

if (nargin < 2)
    if (rem(N,2)==0)
        center = N/2+1;
    else
        center = (N+1)/2;
    end
    P = [1 center 0 0 sqrt(N/4/pi)];
end

if (size(P,2) ~= 5)
    error('Matrix P has the wrong number of columns.')
end

x = zeros(N,1);
% n = (1:N)';
for p = 1:size(P,1)
    A = P(p,1);
    t = P(p,2);
    f = P(p,3);
    cr = P(p,4);
    d = P(p,5);
    clet = chirplet(N,t,f,cr,d,dflag);
    x = x+A*clet;
end

end % function

% -----------------------------
% Subroutines
% -----------------------------
function clet = chirplet(N,t,f,cr,d,dflag)
% CHIRPLET -- Build one unitary chirplet
%
% Sytax:
%   clet = chirplet(N, t,f,cr,d,dflag)
%
% Inputs:
%   N       the signal length
%   t       time center (sample or unit sampling time)
%   f       frequency center (rad)
%   cr      chirp rate (rad/sample)
%   d       time duration of the chirplet (sample or unit sampling time)
%   dflag   consider the discretization effect (Optional: default 0)
%
% Outputs:
%   clet    the unitary chirplet
%

narginchk(5,6);

rep = 5;        % control the accuracy of discretization

if nargin < 5, dflag = 0; end
n = (1:N)';

if dflag        % consider effect
    dcp = zeros(N,1);
    for r = -rep:rep
        acp = exp(-((n+r*N-t)/2/d).^2);
        bcp = exp(1i*cr/2*(n+r*N-t).^2);
        ccp = exp(1i*f*(n+r*N-t));
        dcp = dcp+acp.*bcp.*ccp;
    end %for
    clet = dcp/norm(dcp);
else            % don't care
    am = exp(-((n-t)/2/d).^2) * sqrt(1/sqrt(2*pi)/d); % for normalization 
    chirp = exp(1i * (cr/2*(n-t).^2 + f*(n-t)));
    clet = am.*chirp;
end %if

end % function

% [EOF]
