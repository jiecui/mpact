function [P, res] = mle_adapt_chirplets(x, Q, varargin)
% MLE_ADAPT_CHIRPLETS decompose signal with MLE adaptive chirplet transform
%
%  Syntax:
%       [P, res] = mle_adapt_chirplets(x, Q)
%       [P, res] = mle_adapt_chirplets(x, Q, M)
%       [P, res] = mle_adapt_chirplets(x, Q, M, d)
%       [P, res] = mle_adapt_chirplets(x, Q, M, d, cr)
%       [P, res] = mle_adapt_chirplets(x, Q, M, d, cr, verbose)
%       [P, res] = mle_adapt_chirplets(x, Q, M, d, cr, verbose, mnits)
%       [P, res] = mle_adapt_chirplets(x, Q, M, d, cr, verbose, mnits, level)
%       [P, res] = mle_adapt_chirplets(____, 'RefineAlgorithm', ref_alg)
%
%  Inputs:
%       x       - signal
%       Q       - number of chirplets to look for (if Q = 0, until press
%                 'q' to quit)
%       M       - resolution for Newton-Raphson refinement (optional, 
%                 default = 64)
%       d       - initial value of duration (default = 50)
%       cr      - initial value of chirprate (default = 0)
%       verbose - verbose flag, 'yes', 'no', 'vv' (default = yes)
%       mnits   - maximum number of iterations (optional, default = 5)
%       level   - level of difficulty of MLE
%       ref_alg - { 'expectmax', 'maxlikeliest' }
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
% $ Revision: 0.1 $  $ Date: Thu 02/23/2017  3:48:56.670 PM $
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
d   = p.Results.d;
cr  = p.Results.cr;
verbose = p.Results.verbose;
mnits   = p.Results.mnits;
level   = p.Results.level;
ref_alg = p.Results.RefineAlgorithm;

% =========================================================================
% find chirplets use MP-EM algorithm
% =========================================================================
verbose = lower(verbose);
if strcmp(verbose, 'vv')
    vb = true;
else
    vb = false;
end % if

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
        cprintf('Keywords', '\nEstimated chirplet %d:', i);
        cprintf('Keywords', ...
            '\n|A| = %-6.2f, Tc = %-6.2f, Fc = %-4.2f, Cr = %-7.4f, Dt = %-6.2f\n',...
            abs(P_i(1)), P_i(2), P_i(3), P_i(4), P_i(5));
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
        cprintf('Keywords', '\n--------------');
        cprintf('Keywords', '\n Results');
        cprintf('Keywords', '\n--------------');
        cprintf('Keywords', '\nNo.\t%s\t%s\t%s\t%s\t%s', ...
            '|A|', 'Tc', 'Fc', 'Cr', 'Dt');
        for k = 1:i
            cprintf('Keywords', '\n%-5d\t%-6.2f\t%-7.2f\t%-4.2f\t%-+7.4f\t%-6.2f',...
                k,abs(P(k,1)),P(k,2),P(k,3),P(k,4),P(k,5));
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

end % function

% =========================================================================
% subroutines
% =========================================================================
function p = check_inputs(x, Q, varargin)

validAlgorithm = { 'expectmax', 'maxlikeliest' };
validVerbose = {'yes', 'no', 'vv'};

p = inputParser;

p.addRequired('x', @isnumeric);
p.addRequired('Q', @isnumeric);
p.addOptional('M', 64, @isnumeric);
p.addOptional('d', 50, @isnumeric);
p.addOptional('cr', 0, @isnumeric);
p.addOptional('verbose', 'yes',...
    @(x) any(validatestring(lower(x), validVerbose)));
p.addOptional('mnits', 5, @isnumeric);
p.addOptional('level', 2, @(x) any([0 1 2 3] == x));
p.addParameter('RefineAlgorithm', 'expectmax',...
    @(x) any(validatestring(lower(x), validAlgorithm)));

p.parse(x, Q, varargin{:});

end % funciton

% [EOF]
