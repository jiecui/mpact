function varargout = chirpmpd(varargin)
% CHIRPMPD -- Decomposition using matching pursuit with 
% Gaussian chirplet atoms
%
% Syntax:
%   book = chirpmpd(dict,sig,L,i0,a,verbose)
%
% Inputs:
%   1. DICT    - the dictionary of Gussain chirplet atoms
%   2. SIG     - the signal to be decomposed 
%   3. L       - number of atoms to be collected or extracted
%   4. I0      - the first level to rotate the atom, I0 = 1
%   5. A       - the radix of scale, usually A = 2
%   6. verbose - show details
%
% Outputs:
%   1. book    - code book containing decomposed atoms
% 

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 06-Feb-2004
% $Revision: 1.0$  $Date: 06-Apr-2005 22:40:40$

% --------------------------------------------------
% Main function
% --------------------------------------------------

% Parse inputs
error(nargchk(1,6,nargin));

dict = varargin{1};     % dictionary
sig  = varargin{2};     % the signal need to be analyzed
L    = varargin{3};     % number of atoms
i0   = varargin{4};     % the first level to rotate the atom
a    = varargin{5};     % the radix of scale
if nargin < 6
    verbose = 1;        % default: verbose
else
    verbose = varargin{6};
end %if

% Locat parameters setting
% Initializing book structure
book = struct('code',zeros(L,5),'size',L,'curpoint',1,'corratio',zeros(L,1));
book.N = length(sig);       % set the length of the signal
book.i0 = i0;               % the first level to rotate the atom
book.rad = a;               % the radix of scale

% Error checking

% Main body
% Initialization of the residual energy
Rf0 = sig;              % Rf stands for the residue energy at each iteration
                        % Rf0 is the energy residue at step 0,
                        % initialization, which should be equal to the
                        % signal itself. Note that the signal (sig) can be
                        % a complex signal.
% Initialization of chirpmpd
if verbose, fprintf('\nInitializing, please wait ...'); end %if

[Rf0gbetal,betal,Rf0gbetabar] = localInitMPD(dict,sig,i0,a); 
                                       % gbetal = the atom that has the biggest ampitude; 
                                       % Rf0betabar = <Rf0,gbetabar>, the
                                       % projection of Rf0 on ALL atoms
                                       % gbetabar.
cr = abs(Rf0gbetal)/norm(Rf0);              % cr = correlation coefficients, 
                                            % Rf0gbetal = the projection of
                                            % Rf0 on the atom gbetal
nrmRf1 = sqrt(norm(Rf0)^2-abs(Rf0gbetal)^2);
if verbose
    fprintf('\nNow found the first code, writing into code book ...');
end %if
book = writebook(book,Rf0gbetal,betal,cr);  % write the code beta_l with the highest amplitude 
                                            % into the book

Rf1gbetabar = zeros(size(Rf0gbetabar));     % Rf1gbetabar = Rf1 projected on all atoms gbetabar
                                            % for the next iteration
for l = 1:L-1,
    if verbose
        fprintf('\nNow looking for the %d code ...',l+1);
    end %if
    inprod = atominprod(dict,betal,i0,a);   % inner product of atom betal with all other atoms
    Rf1gbetabar = Rf0gbetabar - Rf0gbetal*inprod;   % UPDATING FORMULA
    % Find the atom with the maximum amplitude
    aR = abs(Rf1gbetabar);                  % amplitudes of projections
    I = find(aR == max(aR(:)));
    [seq,q1,p1] = ind2sub(size(Rf1gbetabar),I);
    Rf1gbetal = Rf1gbetabar(seq(1),q1(1),p1(1));    % find the projection
    [k,m] = seq2idx(seq(1),i0,a);
    betal = [k,m,q1(1)-1,p1(1)-1];
    cr = abs(Rf1gbetal)/nrmRf1;             % correlation coefficients
    book = writebook(book,Rf1gbetal,betal,cr);
    % Update
    nrmRf1 = sqrt(nrmRf1^2-abs(Rf1gbetal)^2); 
    Rf0gbetal = Rf1gbetal;
    Rf0gbetabar = Rf1gbetabar;
end%for

% Outputs
varargout{1} = book;

return;

%--------------------------------------------------------------
function varargout = localInitMPD(varargin)

% Syntax:
%   [Rf0betal,betal,Rf0gbetabar] = localInitMPD(dict,sig,i0,a)

% Parsing the inputs
dict = varargin{1};
Rf0 = varargin{2};
i0 = varargin{3};
a = varargin{4};

% Error checking

% Parameter setting
Rf0 = Rf0(:);              % force to be column vector
nidx = size(dict,2);       % number of atoms in the dictionary
N = length(Rf0);           % the length of the signal
n = 0:N-1; n = n(:);
maxcoe = 0; maxabscoe = 0; % the max coefficient
betal = zeros(1,4);        % the code corresponding to the max
Rf0gbetabar = zeros(nidx,N,N);   % store <Rf0,gbetabar>
Rf0gkm = zeros(N);

% Main boday
% Now matching the atoms one by one, seq: atom sequence (k,m);q: time; p:
% frequency
for seq = 1:nidx,
%     disp(sprintf('trying atom %d ...',seq));
    for q = 0:N-1,
        gkm = dict(:,seq);
        % To do: can use 'repmat'
        tshift = circshift(gkm,q);
        Rf0gkm(:,q+1) = Rf0.*conj(tshift);
    end%for
    Rf0gbeta = fft(Rf0gkm);    % freq by timeshift
    aR = abs(Rf0gbeta);  % p_by_q
    [p1,q1] = find(aR == max(aR(:))); 
    if aR(p1,q1)>maxabscoe
        maxabscoe = aR(p1(1),q1(1));  % if more than one are maximum,get 1st
        maxcoe = Rf0gbeta(p1(1),q1(1));
        [k,m] = seq2idx(seq,i0,a);
        betal = [k,m,q1(1)-1,p1(1)-1];
    end%
    Rf0gbetabar(seq,:,:) = Rf0gbeta.'; % q_by_p, timeshift by frequency
end%for

% Outputs
varargout{1} = maxcoe;
varargout{2} = betal;
varargout{3} = Rf0gbetabar;

return;
