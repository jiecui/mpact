function x = make_chirplets(N, varargin)
% MAKE_CHIRPLETS Construct summation of chirplets
%
% Sytax:
%   x = make_chirplets(N)
%   x = make_chirplets(____, P)
%   x = make_chirplets(____, 'PeriodEff', p_eff)
%   x = make_chirplets(____, 'PType', p_type)
%
% Inputs:
%   N           - required length of signal
%   P           - optional M x 5 matrix of parameters, where M is the number
%                 of chirplets and each row is [A t f cr d] (units depending
%                 on parameter 'Equation'.
%                   A:  complex amplitude A = |A|e^{j\phi}
%                   t:  time center (sample or unit sampling time)
%                   f:  frequency center (rad)
%                   cr: chirprate (rad/sample)
%                   d:  chirplet duration (sample)
%                   (default [1, N/2, 0, 0, sqrt(N/4/pi)])
%   PeriodEff   - parameter of discretization effect, p_eff = {true, false}
%                 (default p_eff = true)
%   PType       - parameter of the type of P, p_type = {'Oneill', 'Cohen'}
%                 (default 'ONeill')
% 
% Outputs:
%   x           - constructed signal
%
% Note:
%   In O'Neill's equation, 'd' is the standard deviation of the guassian. d
%   = \sqrt{\frac{N}{4\pi}} gives an atom with a circular Wigner
%   distribution, and 2*sqrt(2)*d is the Rayleigh limit. Assume sampling
%   frequency = 1 Hz and Nyquist = \pi. Use rad for frequency unit.
%
%   In Cohen's equation, the chirplet formula is A(\frac{d}{\pi})^{1/4}
%   e^{-\frac{d}{2}(n-t)^2} e^{j[\frac{cr}{2}(n-t)-f](n-t)}. Assume
%   sampling frequency is N and Nyquist = \pi.
% 
% Examples:
%   N = 100; x = make_chirplets(N, [1 N/2 pi/2 pi/N sqrt(2*N/pi)]); 
%   N = 100; x = make_chirplets(N, [1, 1/2, pi/2, pi, pi/4], 'PType', 'Cohen');
% 
% References:
%   [1] Cohen, L. (1995). Time-frequency analysis. Englewood Cliffs, N.J,
%       Prentice Hall PTR.
%   [2] O'Neill, J. C. and P. Flandrin (1998). Chirp hunting. Proceedings 
%       of the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis (Cat. No.98TH8380).
% 
% See also .

% Copyright 2005-2017 Richard J. Cui. Created: Thu 03/02/2005 10:33:22.814 AM;
% $ Revision: 0.9 Date: Mon 05/29/2017  4:14:43.499 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% -------------------------------------------------------------------------
% parse the inputs
% -------------------------------------------------------------------------
p = check_inputs(N, varargin{:});

N       = p.N;
P       = p.P;
p_eff   = p.p_eff;
p_type  = p.p_type;

% -------------------------------------------------------------------------
% convert chirplet parameters to O'Neill's if necessary
% -------------------------------------------------------------------------
P = convert_p(N, P, p_type);

% -------------------------------------------------------------------------
% construct chirplets
% -------------------------------------------------------------------------
% use O'Neill's equation
x = zeros(N, 1);
for k = 1:size(P,1)
    A_k = P(k, 1);
    t_k = P(k, 2);
    f_k = P(k, 3);
    cr_k = P(k, 4);
    d_k = P(k,5);
    clet_k = chirplet(N, t_k, f_k, cr_k, d_k, p_eff);
    x = x+A_k*clet_k;
end

end % function

% =========================================================================
% subroutines
% =========================================================================
function P = convert_p(N, P_in, p_type)
% convert Cohen paras to O'Neill

switch p_type
    case 'oneill'
        P = P_in;
    case 'cohen' % --> O'Neill
        P = PCo2On(N, P_in);
    otherwise
        error('make_chirplets:convert_p', 'Unknown P type.')
end % switch

end % function

function y = isValidChirplets(P)
% check if valid chirplet parameters

if size(P,2) ~= 5
    y = false;
else
    y = true;
end

end % funciton

function q = check_inputs(N, varargin)
% parse inputs
valid_eq = {'oneill', 'cohen'};

p = inputParser;

p.addRequired('N', @isnumeric);
p.addOptional('P', [1 N/2 0 0 sqrt(N/4/pi)], @isValidChirplets); % check if valid chirplet parameters
p.addParameter('PeriodEff', true, @islogical);
p.addParameter('PType', 'ONeill', @(x) any(validatestring(lower(x), valid_eq)));

p.parse(N, varargin{:});

q.N     = p.Results.N;
q.P     = p.Results.P;
q.p_eff = p.Results.PeriodEff;
q.p_type= lower(p.Results.PType);

end % function

function clet = chirplet(N,t,f,cr,d,dflag)
% chirplet Build one unitary chirplet
%
% Sytax
%   clet = chirplet(N, t,f,cr,d,dflag)
%
% Inputs
%   N       - the signal length
%   t       - time center (sample or unit sampling time)
%   f       - frequency center (rad)
%   cr      - chirp rate (rad/sample)
%   d       - time duration of the chirplet (sample or unit sampling time)
%   dflag   - consider the discretization effect (Optional: default 0)
%
% Outputs
%
%   clet    the unitary chirplet

narginchk(5,6);

rep = 5;        % control the accuracy of discretization

if nargin < 5, dflag = 0; end
n = (1:N)';

if dflag == true % consider effect
    dcp = zeros(N,1);
    for r = -rep:rep
        acp = exp(-((n+r*N-t)/2/d).^2);
        bcp = exp(1i*cr/2*(n+r*N-t).^2);
        ccp = exp(1i*f*(n+r*N-t));
        dcp = dcp+acp.*bcp.*ccp;
    end %for
    clet = dcp/norm(dcp); % ?
else % don't consider periodization
    am = exp(-((n-t)/2/d).^2) * sqrt(1/sqrt(2*pi)/d); % for normalization
    chirp = exp(1i * (cr/2*(n-t).^2 + f*(n-t)));
    clet = am.*chirp;
end %if

end % function

% [EOF]