function book = max_chirpmpd(x, D, i0, a, verbose)
% MAX_CHIRPMPD find the maximumn chirplet using matching pursuit with Gaussian chirplet atoms
%
% Syntax:
%   book = max_chirpmpd(x, D, i0, a)
%   book = max_chirpmpd(x, D, i0, a, verbose)
%
% Inputs:
%   1. x        - the signal to be decomposed
%   2. D        - the depth of the dictionary
%   3. i0       - the first level to rotate the atom, I0 = 1
%   4. a        - the radix of scale, usually A = 1
%   5. verbose  - show details (default: verbose = 1)
%
% Outputs:
%   1. book    - code book containing decomposed atoms
%

% Copyright 2005-2016 Richard J. Cui. Created: Fri 04/08/2005 10:33:22.814 AM
% $ Revision: 0.5 $  $ Date: Wed 12/14/2016  2:59:45.038 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com


% =========================================================================
% Input parameters and options
% =========================================================================
narginchk(1,5);

if nargin < 5
    verbose = true; % default
end %if

% =========================================================================
% matching-pursuit
% =========================================================================
% Initializing book structure
% ---------------------------
L = 1;
book = struct('code',zeros(L,5),'size',L,'curpoint',1,'corratio',zeros(L,1));
book.N  = length(x); % set the length of the signal
book.i0 = i0; % the first level to rotate the atom
book.rad = a; % the radix of scale

% Error checking

% Initialization of the residual energy
% -------------------------------------
% Rf0 is the energy residue at step 0, initialization, which should be
% equal to the signal itself. Note that the signal (x) can be a complex
% signal.
Rf0 = x; % Rf stands for the residue energy at each iteration

% --------------
% MP of chirpmpd
% --------------
if verbose
    fprintf('\nLooking for the chirplet with the maximum amplitude...'); 
end %if

% =============================================================
% Rf0gbetal = the projection of Rf0 on the atom gbetal
[Rf0gbetal, betal] = mp_chirplet_bultan(Rf0, D, i0, a, verbose);
% ==============================================================

cc = abs(Rf0gbetal)/norm(Rf0);  % cc = correlation coefficients,
if verbose
    fprintf('\nNow found the code, writing into code book...\n');
end %if

% write the code beta_l with the highest amplitude into the book
book = writebook(book, Rf0gbetal, betal, cc);

end % function

% =========================================================================
% subroutines
% =========================================================================

% [EOF]