function P = refinechirp(P0, vlb, vub, x, verbose)
% REFINECHIRP Refine the chirplet using Qausi-Newton method
% 
%  Syntax:
%       P = refinechirp(theta0, vlb, vub, x, verbose)
%
%  Inputs:
%       P0      - intitial estimate of chirplet parameters (from Matching
%                 Pursuit results)
%       vlb     - lower bound of the variables
%       vub     - upper bound of the variables
%       x       - the signal to be fitted
%       verbose - whether verbose displays
%
%  Outputs:
%       theta   - vector of refined chirplet parameters (see make_chirplets.m)
%

% Copyright 2005-2016 Richard J. Cui. Created: Thr 02/22/2005 10:33:22.814 AM
% $ Revision: 0.2 $  $ Date: Tue 12/06/2016  6:01:18.628 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com


% =========================================================================
% Set options
% =========================================================================
opt = optimset;
opt.Display = 'off';        % Display no output.
opt.Diagnostics = 'off';    % Print diagnoistic information
opt.GradObj = 'on';         % gradient of objective is provided
opt.MaxFunEvals = 1000;     % Maximum number of function evaluations allowed
opt.MaxIter = 1000;         % Maximum number of iterations allowed
opt.TolFun = 1e-6;          % Termination tolerance on the function value.
opt.TolX = 1e-6;            % Termination tolerance on x.
if verbose
    opt.Display = 'iter';   % Show details.
end%if

% =========================================================================
% optimize by Newton method
% =========================================================================
f = @(P)obj_chirp(P, x);
P = fmincon(f, P0, [], [], [], [], vlb, vub, [], opt);

end

% [EOF]