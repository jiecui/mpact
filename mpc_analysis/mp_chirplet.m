function P = mp_chirplet(x, varargin)
% MP_CHIRPLET Estimate a chirplet that best fits the signal using Matching Pursuit algorithm
%
%  Syntax:
%   P = mp_chirplet(x)
%   P = mp_chirplet(x, M)
%   P = mp_chirplet(x, M, D)
%   P = mp_chirplet(x, M, D, i0)
%   P = mp_chirplet(x, M, D, i0, radix)
%   P = mp_chirplet(x, M, D, i0, radix, verbose)
%
%  Inputs:
%   x           - signal
%   M           - resolution for Newton-Raphson refinement (optional, 
%                 default is 64)
%   D           - The depth of decomposition
%   i0          - The first level to rotate the chirplets
%   radix       - radix of scale
%   verbose     - verbose flag (optional, default = false)
%
%  Outputs:
%   P           - vector of chirplet parameters (see make_chirplets.m)
% 
% Note:
%   The estimated chirplet is locally refined with Newton-Raphson method.
%   Assume sampling frequency of signal is 1 Hz.
% 
% Example:
%
% References:
%   [1]	S. G. Mallat and Z. Zhang, "Matching pursuit with time-frequency
%       dictionaries," IEEE Transactions on Signal Processing, vol. 41,
%       pp. 3397-3415, December 1993.
%   [2]	A. Bultan, "A four-parameter atomic decomposition of chirplets,"
%       IEEE Transactions on Signal Processing, vol. 47, pp. 731-745, Mar
%       1999.
%   [3]	J. C. O'Neill and P. Flandrin, "Chirp hunting," in Proceedings of 
%       the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis, 1998, pp. 425-428.
%
% See also make_chirplets.

% Copyright 2005-2016 Richard J. Cui. Created: Tue 02/22/2005 10:33:22.814 AM
% $ Revision: 1.1 $  $ Date: Wed 12/14/2016 11:30:38.176 AM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% =========================================================================
% Input parameters and options
% =========================================================================
p = check_inputs(x, varargin{:});
x   = p.Results.x;
M   = p.Results.M;
D   = p.Results.D;
i0  = p.Results.i0;
radix = p.Results.radix;
verbose = p.Results.verbose;

% this works better numerically
x = x(:); % force to be column vector
e = norm(x);
x = x/e;
N = length(x);

% =========================================================================
% estimate the chirplet globally: MP + Newton-Raphson
% =========================================================================
% ------------------------------
% estimate the chirplet using MP
% ------------------------------
if verbose
    fprintf('\n---------------------------');
    fprintf('\nEstimates of initial values');
    fprintf('\n---------------------------');
end %if

% MPing the chirplet...
new_book = max_chirpmpd(x, D, i0, radix, verbose);
[tc, fc, cr, dt] = bultan2paras(new_book, N, i0, radix);

% show initial guess if verbose
if verbose
    A = new_book.code(1);
    fprintf('\nInitial estimates:\n|A| = %-6.2f, Tc = %-6.2f, Fc = %-4.2f, Cr = %-7.4f, Dt = %-6.2f\n',...
        abs(A), tc, fc, cr, dt);
end %if

% ----------------------------------------------
% Do a quasi-Newton optimization on the chirplet
% ----------------------------------------------
if verbose
    fprintf('\n------------------------------------');
    fprintf('\n Newton-Raphson refinement (FMINCON)');
    fprintf('\n------------------------------------');
    fprintf('\n')
end %if

% newton-raphson refinement
P = nr_refine(x, e, M, tc, fc, cr, dt, verbose);

if verbose
    fprintf('\nRefined estimates:\n|A| = %-6.2f, Tc = %-6.2f, Fc = %-4.2f, Cr = %-7.4f, Dt = %-6.2f\n',...
        abs(P(1)), P(2), P(3), P(4), P(5));
end % if

end

% =========================================================================
% subroutines
% =========================================================================
function P = nr_refine(x, e, M, tc, fc, cr, dt, verbose)
% refine the chirplet with Newton-Raphson method

N = length(x);

% a longer window is useful here.
Z = 4;
rt = round(tc);
if (rt - Z * M < 1) && (rt + Z * M > N)
    xx = [zeros(Z * M - rt + 1, 1); x; zeros(Z * M - N + rt, 1)];
elseif rt - Z * M < 1
    xx = [zeros(Z * M - rt + 1, 1); x(1:rt + Z * M)];
elseif rt + Z * M > N
    xx = [x(rt - Z * M:N); zeros(Z * M - N + rt, 1)];
else
    xx = x(rt - Z * M:rt + Z * M);
end

x0 = [Z*M+1+(tc-rt) fc cr dt];
vlb = [1 0 -inf .25];
vub = [2*Z*M+1 2*pi inf N/2];
Pr = refinechirp(x0, vlb, vub, xx, verbose);
Pr(1) = rt + Pr(1) - (Z*M+1);
Pr(2) = mod(Pr(2),2*pi);

y = make_chirplets(N, [1 Pr]);
A = y'*x;
A = A*e;

P = [A Pr];

end % function

function p = check_inputs(x, varargin)

p = inputParser;

p.addRequired('x', @isnumeric);
p.addOptional('M', 64, @isnumeric);
p.addOptional('D', 5, @isnumeric);
p.addOptional('i0', 1, @isnumeric);
p.addOptional('radix', 2, @isnumeric);
p.addOptional('verbose', false, @islogical);

p.parse(x, varargin{:});

end % funciton

function [tc, fc, cr, dt] = bultan2paras(cbook, N, i0, radix)
% obtain chirplet parameters from Bultan indexes

k_idx = cbook.code(2);        % scale index
m_idx = cbook.code(3);        % rotation index
q_idx = cbook.code(4);        % time-shift index
p_idx = cbook.code(5);        % freq-shift index

% refer to Bultan,99 (10) & (11) and O'Neil,98 (0) for Cr and Dt
angm = getalpham(k_idx, m_idx, i0, radix); % get the discrete angle = alpha_m
s = radix^k_idx; % scale_k
sigma = sqrt(sin(angm)^2 ...            % sigma(scale_k, ang_m)
    +s^4*cos(angm)^2)/s;

% -----------------------------------------------------
% time center = Tc, index --> second (assmue fs = 1 Hz)
% -----------------------------------------------------
tc = q_idx;

% -----------------------------------------------------
% freq center = Fc, index --> rad
% -----------------------------------------------------
fc = p_idx / N * 2 * pi;

% index --> chirp rate (Cr) and Duration (Dt)
% -----------------------------------------------------
% chirp rate = Cr
% -----------------------------------------------------
% xi = ((s^4-1)*cos(angm)*sin(angm))...
%     /(sin(angm)^2+s^4*cos(angm)^2);
% c = xi; % xi&sigma Bultan's definitions
cr = 2 * pi / N * tan(angm);
% -----------------------------------------------------
% time spread / duration = Dt
% -----------------------------------------------------
% d = sigma/sqrt(2); % c&d O'Neil's definition
dt = sigma*sqrt(N/pi/2);

end % function

% [EOF]