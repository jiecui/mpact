function P = best_chirplet(x, level, M, t, f, c, d, verbose)
% BEST_CHIRPLET Find the best chirplet with maximum likelihood estimation 
%
%  Syntax:
%   P = best_chirplet(x, level, M, c, d, t, f, verbose)
%
%  Inputs:
%   x           - signal
%   level       - level of difficulty, 0 (easiest) -> 3 (hardest) (optional,
%                 default is 2)
%                   0: quasi-Newton (QN)
%                   1: est_tf -> est_c -> est_d -> QN
%                   2: est_cd_global -> est_tf -> est_c -> est_d -> QN
%                   3: est_cd_global -> (est_tf -> est_c -> est_d) x 3 -> QN
%   M           - resolution parameter (optional, default is 64)
%   c,d,t,f     - possible intializations; level 0 uses c, d, t, and f;
%                 level 1 uses c and d; (optional, defaults are 0, 50, N/2
%                 and 0, respectively)
%   verbose     - verbose flag (optional, default is 0)
%
%  Outputs:
%   P           - vector of chirplet parameters (see make_chirplets.m)
%
%  References
%   [1]	J. C. O'Neill and P. Flandrin, "Chirp hunting," in Proceedings of 
%       the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis, 1998, pp. 425-428.
% 
% Example:
% 
% See also make_chirplets.

% Copyright 2005-2016 Richard J. Cui. Created: Tue 02/22/2005  1:23:52.278 PM
% $ Revision: 0.7 $  $ Date: Thu 12/15/2016 10:18:53.485 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

narginchk(1,8);

x = x(:);
N = length(x);

% this works better numerically
e = norm(x);
x = x/e;

if (nargin<2), level = 2; end
if (nargin<3), M = 64; end;
if (nargin<4), verbose = 0; end;
if (nargin<5), c = 0; end;
if (nargin<6), d = 50; end;
if (nargin<7), t = N/2; end;
if (nargin<8), f = 0; end;

if any([0 1 2 3] == level) == false
    error('choose a valid level')
end

% estimate chirp-rate and duration globally
if verbose
    fprintf('\n---------------------------');
    fprintf('\nEstimates of initial values');
    fprintf('\n---------------------------');
end %if
r = 5; % robustness parameter
if (level == 2 || level == 3)
    [c, d] = est_cd_global(x, M, r);
    if verbose
        fprintf('\n');
        fprintf(1,'est_cd_global -> c = %7.4f, d = %6.2f\n', c, d);
    end
end

if (level == 1 || level == 2)
    [t, f] = est_tf(x, c, d, M);
    
    f = mod(f, 2*pi);
    if verbose, fprintf(1,'est_tf -> t = %7.2f, f = %4.2f\n', t, f); end
    
    c = est_c(x, t, f, M);
    if verbose, fprintf(1,'est_c -> c = %7.4f\n', c); end
    
    d = est_d(x, t, f, c, M);
    if verbose, fprintf(1,'est_d -> d = %6.2f\n', d); end
elseif (level == 3)
    for i=1:3
        [t, f] = est_tf(x, c, d, M);
        f = mod(f, 2*pi);
        %if verbose, fprintf(1,'est_tf -> t = %7.2f, f = %4.2f\n', t, f); end
        c = est_c(x, t, f, M);
        %if verbose, fprintf(1,'est_c -> c = %7.4f\n', c); end
        d = est_d(x, t, f, c, M);
        %if verbose, fprintf(1,'est_d -> d = %6.2f\n', d); end
    end
end
% show initial guess if verbose
if verbose
    fprintf('\nTc = %-6.2f, Fc = %-4.2f, Cr = %-+7.4f, Dt = %-6.2f\n',...
        t,f,c,d);
end %if

% Do a quasi-Newton maximization on the windowed signal
if verbose
    fprintf('\n-----------------------------------');
    fprintf('\n Quasi-Newton refinement (FMINCON)');
    fprintf('\n-----------------------------------');
    fprintf('\n');
end %if
% a longer window is useful here.
Z = 4;
rt = round(t);
if ( (rt-Z*M < 1) && (rt+Z*M > N) )
    xx = [zeros(Z*M-rt+1,1) ; x ; zeros(Z*M-N+rt,1)];
elseif (rt-Z*M < 1)
    xx = [zeros(Z*M-rt+1,1) ; x(1:rt+Z*M)];
elseif (rt+Z*M > N)
    xx = [x(rt-Z*M:N) ; zeros(Z*M-N+rt,1)];
else
    xx = x(rt-Z*M:rt+Z*M);
end

x0 = [Z*M+1+(t-rt) f c d];
vlb = [1 0 -inf .25];
vub = [2*Z*M+1 2*pi inf N/2];
P = refinechirp(x0,vlb,vub,xx,verbose);
P(1) = rt + P(1) - (Z*M+1);
P(2) = mod(P(2),2*pi);
if verbose
    fprintf('\nTc = %-7.2f, Fc = %-4.2f, Cr = %-+7.4f, Dt = %-6.2f\n',...
        P(1),P(2),P(3),P(4));
end %if

y = make_chirplets(N, [1 P]);
A = y'*x;
A = A*e;

P = [A P];

end % function

% [EOF]