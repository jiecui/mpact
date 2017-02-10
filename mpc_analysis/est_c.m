function c = est_c(x,t,f,M,method)
% EST_C -- Estimate the chirp rate from a local measure.
%
%  Syntax
%    c = est_c(x,t,f,M,method)
%
%  Inputs
%    x       signal vector
%    t       current estimate of time location
%    f       current estimate of frequency location
%    M       resolution parameter (optional, default is 64)
%    method  use the Wigner distribution 'wig' or local ambiguity
%            function 'laf' (optional, default is 'wig')
%
%  Outputs
%    c      estimate of the chirp rate
%
% The 'laf' method is not yet implemented.

% Author(s):    JOC, R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; 02-Mar-2005
% $Revision: 0.1 $  $Date: 02-Mar-2005 22:25:16$

narginchk(3, 5);

if (nargin < 4)
  M = 64;
end
if (nargin < 5)
  method = 'wig';
end

M = 4*floor(M/4) ;  % want M to be a multiple of 4
x = x(:);
N = length(x);

% center in time and window, window length = M
% should probably use a non-rectangular window -> gaussian, hamming etc.
MM = M/2;
rt = round(N/2);
if M > N
    xx = zeros(M, 1);
    t0 = MM - rt;
    t1 = t0 + N -1;
    xx(t0:t1) = x;
    x = xx;
elseif M < N
    t0 = rt - MM;
    t1 = t0 + M - 1;
    x = x(t0:t1);
end % if

% center in frequency
x = x .* exp(-1i*(1:M)'*f);

X = fftshift(fft(fftshift(x)))/sqrt(M);

% estimate the chirp rate
R1 = zeros(M/2,1);
R2 = zeros(M/2,1);
nn = (-M/2:1:M/2-1)';
for a = 1:MM
  angle = a*90/MM - 45;
  c = angle2cr(angle,M,2*pi);

  y = x .* exp(-1i*c/2*nn.^2);
  R1(a) = abs(sum(y))^2;  % = \int W_y(t,0) dt

  Y = X .* exp(-1i*c/2*nn.^2);
  R2(a) = M/2/pi*abs(c)*abs(sum(Y))^2;  % = \int W_y(0,omega) d omega
end

R = [R2(M/4+1:end) ; R1 ; R2(1:M/4-1)];

[m, a] = max(R);
angle = a*180/M - 90;
c = angle2cr(angle, M, 2*pi);
