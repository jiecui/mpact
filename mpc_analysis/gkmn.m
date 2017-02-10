function g = gkmn(k, M, i0, a, N)
% GKMN Get Gaussian chirplet atom at certian scale and rotational angle
% 
% Syntax:
%   g = gkmn(N, k, M, i0, a)
%
% Input(s):
%   k   - the index of scale
%   M   - the index of rotation
%   i0  - the first level to rotate the atom, usually I = 1
%   a   - radix of scale, usually a = 2
%   N   - imput signal length, C^N (N dimensional complex number)
%
% Output(s):
%   g   - the Guassian chirplet atom
% 
% Example:
%
% See also mp_chiprlet_bultan.

% Copyright 2005-2016 Richard J. Cui. Created: Thu 12/21/2005 10:33:22.814 AM
% $ Revision: 0.2 $  $ Date: Tue 12/06/2016  6:01:18.628 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% --------------------------------------------------
% Main function
% --------------------------------------------------

% Main body
ang_m = getalpham(k, M, i0, a); % get the discrete angle = alpha_m
csig  = sigkm(a^k, ang_m);         % current sigma(scale_k, ang_m)
cx    = xkm(a^k, ang_m);           % current x(scale_k,ang_m)
cg    = gkm(csig, cx, N);          % g(sig_km,xi_km,n)
g     = cg(:); % force to be column

end % function

% =========================================================================
% subroutines
% =========================================================================
function varargout = sigkm(varargin)

% csig = sigkm(a^k,angm)
% sigma

s = varargin{1};
angm = varargin{2};

sig = sqrt(sin(angm)^2+s^4*cos(angm)^2)/s;
varargout{1} = sig;

end % function

%------------------------------------------------------------
function varargout = xkm(varargin)

% cx = xkm(a^k,angm)
% xi

s = varargin{1};
angm = varargin{2};

x = ((s^4-1)*cos(angm)*sin(angm))/(sin(angm)^2+s^4*cos(angm)^2);
varargout{1} = x;

end % function

%------------------------------------------------------------
function varargout = gkm(varargin)

% cg = gkm(csig,cx,N)
% Gussain atom

csig = varargin{1};
cx = varargin{2};
N = varargin{3};

d = 5;          % this control the accuracy of the Gaussian window

cg = zeros(1,N);

% assume delta_t = sqrt(2*pi/N)

% Consider the effect of digitizing
for n = 0:N-1
    for r = -d:d
        cg(n+1) = cg(n+1)+exp(-pi/N*(1/csig^2-1i*cx)*(n+r*N)^2);
    end%for
end%for

% No consideration of digitizing
% cg = sqrt(1/sqrt(pi)*csig)*exp(-(1/csig^2-j*cx)*(sqrt(2*pi/N)*(0:N-1)).^2/2);

% Normailization
cg = cg/norm(cg);

varargout{1} = cg;

end % function

% [EOF]