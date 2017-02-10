function x = fmfun(N,A,p,phs)
% fmfun -- create a signal with a sinusoidal frequency modulation
%
%  Usage
%    x = fmsin(N,A,p,phs)
%
%  Inputs
%    N    length of the signal
%    A    amplitude of frequency modulation
%    p    number of cycles (optional, default is 1)
%    phs  phase shift of the instantaneous frequency (optional, default is 0)
%
%  Outputs
%    x    signal

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 13-Oct-2005
% $Revision: 0.1 $  $Date: 13-Oct-2005 14:46:05$

error(nargchk(1,4,nargin));

if (nargin < 4)
    phs = 0;
end
if (nargin < 3)
    p = 1;
end
if (nargin < 2)
    A = 1;
end

n = (0:N-1)';
x = exp(-j*N/4/p*A*cos(p*2*pi*n/N + phs));
