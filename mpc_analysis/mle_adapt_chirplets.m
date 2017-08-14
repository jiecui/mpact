function [P, res] = mle_adapt_chirplets(x, Q, varargin)
% MLE_ADAPT_CHIRPLETS decompose signal with MLE adaptive chirplet transform
%
%  Syntax:
%       [P, res] = mle_adapt_chirplets(x, Q)
%       [P, res] = mle_adapt_chirplets(____, M)
%       [P, res] = mle_adapt_chirplets(____, M, CP0)
%       [P, res] = mle_adapt_chirplets(____, M, CP0, verbose)
%       [P, res] = mle_adapt_chirplets(____, M, CP0, verbose, mnits)
%       [P, res] = mle_adapt_chirplets(____, M, CP0, verbose, mnits, level)
%       [P, res] = mle_adapt_chirplets(____, 'RefineAlgorithm', ref_alg)
%       [P, res] = mle_adapt_chirplets(____, 'PType', p_type)
%
%  Inputs:
%       x       - signal
%       Q       - number of chirplets to look for (if Q = 0, until press
%                 'q' to quit)
%       M       - resolution for Newton-Raphson refinement (optional, 
%                 default = 64)
%       CP0     - initial value of chirplet [tc, fc, cr, d] (default = [1 pi/2 0 1])
%       verbose - verbose flag, 'yes', 'no', 'vv' (default = yes)
%       mnits   - maximum number of iterations (optional, default = 5)
%       level   - level of difficulty of MLE
%       ref_alg - { 'expectmax', 'maxlikeliest' }
%       PType   - (parameter) type of chirplet parameter P, p_type = 
%                 {'Oneill', 'Cohen'}
%
%  Outputs:
%       P       - Q_by_5 matrix of chirplet parameters (see make_chirplets.m)
%       res     - norm of the signal and the residuals for 1 to Q chirplets;
%                 could be used for a selection of stopping criterion
%
% Example:
%
% References:
%   [1]	J. Cui and W. Wong, "The adaptive chirplet transform and visual 
%       evoked potentials," IEEE Transactions on Biomedical Engineering,
%       vol. 53, pp. 1378-1384, Jul 2006.
%   [2]	J. C. O'Neill and P. Flandrin, "Chirp hunting," in Proceedings of
%       the IEEE-SP International Symposium on Time-Frequency and
%       Time-Scale Analysis, 1998, pp. 425-428.
% 
% See also make_chirplets.

% Copyright 2017 Richard J. Cui. Created: Thu 02/23/2017  3:48:56.670 PM
% $ Revision: 0.2 $  $ Date: Wed 05/31/2017  4:17:25.467 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% =========================================================================
% parse inputs
% =========================================================================
p = check_inputs(x, Q, varargin{:});

x   = p.Results.x;
Q   = p.Results.Q;
M   = p.Results.M;
CP0   = p.Results.CP0;
verbose = p.Results.verbose;
mnits   = p.Results.mnits;
level   = p.Results.level;
ref_alg = p.Results.RefineAlgorithm;
p_type  = lower(p.Results.PType);

N = length(x);

% =========================================================================
% find chirplets use MP-EM algorithm
% =========================================================================
% verbose flag
verbose = lower(verbose);
if strcmp(verbose, 'vv')
    vb = true;
else
    vb = false;
end % if

% initial parameters
tc  = CP0(1);
fc  = CP0(2);
cr  = CP0(3);
d   = CP0(4);

x = x(:);
res = norm(x);

P = [];
done = false;
i = 1;
e = x; % residual
while done == false
    % --------------------------------------------
    % find a new chirplet with MLE
    % --------------------------------------------
    if strcmp(verbose, 'yes') || strcmp(verbose, 'vv')
        cprintf('Keywords', '\n******************************************');
        cprintf('Keywords', '\nSingle chirplet estimation - MLE algorithm');
        cprintf('Keywords', '\n******************************************');
    end % if
    
    % esimtate the i-th chirplet
    P_i = best_chirplet(e, level, M, tc, fc, cr, d, vb);
    P = cat(1, P, P_i);
    
    if strcmp(verbose, 'yes') || strcmp(verbose, 'vv')
        % display P in correct format
        switch p_type
            case 'cohen'
                P_disp = POn2Co(N, P_i);
            case 'oneill'
                P_disp = P_i;
        end % switch
        cprintf('Keywords', '\nEstimated chirplet %d:', i);
        cprintf('Keywords', ...
            '\n|A| = %-6.2f, Tc = %-6.2f, Fc = %-4.2f, Cr = %-7.4f, D = %-6.2f\n',...
            abs(P_disp(1)), P_disp(2), P_disp(3), P_disp(4), P_disp(5));
    end % if
    % ------------------------------
    % refine the estimated chirplets 
    % ------------------------------
    switch lower(ref_alg)
        case 'expectmax' % Expectation-Maximization (EM)
            if (strcmp(verbose, 'yes') || strcmp(verbose, 'vv')) && i > 1
                % show roadmaps
                cprintf('Keywords', '\n********************************************');
                cprintf('Keywords', '\nMultiple chirplets refinement - EM algorithm');
                cprintf('Keywords', '\n********************************************');
                if strcmp(verbose, 'vv')
                    fprintf(sprintf('\nEM iterations(max = %d): ', mnits)); 
                end
            end %if
            [P, e, res] = em_chirplets(x, P, res, M, D, i0, radix, mnits, vb); % EM
        case 'maxlikeliest' % Maximum-liklihood estimation (MLE)
            if (strcmp(verbose, 'yes') || strcmp(verbose, 'vv')) && i > 1
                % show roadmaps
                cprintf('Keywords', '\n*********************************************');
                cprintf('Keywords', '\nMultiple chirplets refinement - MLE algorithm');
                cprintf('Keywords', '\n*********************************************');
                if strcmp(verbose, 'vv')
                    fprintf(sprintf('\nMLE iterations(max = %d): ', mnits)); 
                end
            end %if
            [P, e, res] = mle_chirplets(x, P, res, level, M, mnits, vb); % MLE
    end % switch
    
    % show details if asked for...
    if (strcmp(verbose, 'yes') || strcmp(verbose, 'vv')) && i > 1
        % display P in correct format
        switch p_type
            case 'cohen'
                P_disp = POn2Co(N, P);
            case 'oneill'
                P_disp = P;
        end % switch
        cprintf('Keywords', '\n--------------');
        cprintf('Keywords', '\n Results');
        cprintf('Keywords', '\n--------------');
        cprintf('Keywords', '\nNo.\t%s\t%s\t%s\t%s\t%s', ...
            '|A|', 'Tc', 'Fc', 'Cr', 'D');
        for k = 1:i
            cprintf('Keywords', '\n%-5d\t%-6.2f\t%-7.2f\t%-4.2f\t%-+7.4f\t%-6.2f',...
                k,abs(P_disp(k,1)),P_disp(k,2),P_disp(k,3),P_disp(k,4),P_disp(k,5));
        end
        fprintf('\n');
    end %if
    
    % termination criteria - press 'q' to quit
    if Q == 0
        [~, ~, button] = ginput(1);     % no use of mouse position at present
        if button == 'q', done = true; end
    elseif Q == i
        done = true;
    end %if
    i = i + 1;
end % while

% convert P format if necessary
switch p_type
    case 'cohen'
        P = POn2Co(N, P);
end % switch

end % function

% =========================================================================
% subroutines
% =========================================================================
function yn = check_cp0(x)

yn = isnumeric(x) && length(x) == 4;

end % function

function p = check_inputs(x, Q, varargin)

validAlgorithm = { 'expectmax', 'maxlikeliest' };
validVerbose = {'yes', 'no', 'vv'};
ptype_str = {'oneill', 'cohen'};

p = inputParser;

p.addRequired('x', @isnumeric);
p.addRequired('Q', @isnumeric);
p.addOptional('M', 64, @isnumeric);
p.addOptional('CP0', [1 pi/2 0 1], @check_cp0);
p.addOptional('verbose', 'yes',...
    @(x) any(validatestring(lower(x), validVerbose)));
p.addOptional('mnits', 5, @isnumeric);
p.addOptional('level', 2, @(x) any([0 1 2 3] == x));
p.addParameter('RefineAlgorithm', 'expectmax',...
    @(x) any(validatestring(lower(x), validAlgorithm)));
p.addParameter('PType', 'Oneill', @(x) any(validatestring(lower(x), ptype_str)));

p.parse(x, Q, varargin{:});

end % funciton

% [EOF]
