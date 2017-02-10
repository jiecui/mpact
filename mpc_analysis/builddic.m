function varargout = builddic(varargin)

% BUILDDIC Build dictionary of Gaussian chirplet atoms
% DIC = BUILDDIC(N,A,D,I0)
% Input parameters:
% 1. N        - the dimension of the input signal, beloning to the set of C^N
% 2. A        - the radix of scale, usually A = 2
% 3. D        - the depth of decomposition, scale level
% 4. I0       - the first level to rotate the atom, usually I0 = 1
% Output parameters:
% 1. DIC      - the dictionary containing atoms of Guassian chirplets, N x Nidx (total # atoms)
% 
% Update: Jun. 27, 2004
% Update: Jan. 4, 5, 8, 2004
% 
% Richard J. Cui
% Sensory communication lab
% IBBME, U of Toronto
% 4 Taddle Creek Road
% Toronto, ON M5S 3G9
% Canada
% Email: richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 4-Jan-2004
% $Revision: 0.5 $  $Date: 19-Oct-2005 16:32:26$


% --------------------------------------------------
% Main function
% --------------------------------------------------

% Parse inputs
N   = varargin{1};          % the size of the input signal
a   = varargin{2};          % the radix of scale
D   = varargin{3};          % the depth of decomposition
i0  = varargin{4};          % the first level to rotate the atom

% Create index and dictionary
r = 0:(D-i0-1);                 % total # of levels, r: 0,1,...,D-i0-1
nidx = i0+sum(4*a.^(2*r)-1);    % total number of atoms, from r

atombook = zeros(N,nidx);
% Calculate the scaling-rotation gkm(n)
% k   - discrete index of scale
% m   - discrete index of rotation
% n   - time position
disp(sprintf('%d atoms will be added into the dictionary ...',nidx));
prog = 0;   % intialize the progress index
% hm = floor(nidx/10);
for seq = 1:nidx,
    prog = showprog(seq,nidx,prog); % showing the progress
    % if mod(seq,hm) == 0, fprintf(1,'* '); end
    [k,m] = seq2idx(seq,i0,a);      % get index of scale and rotation from the sequence
    cg = gkmn(k,m,i0,a,N);          % gkmn: Gaussian chirplet atom at scale k and rotation m
    atombook(:,seq) = cg;           % write it into atom book
end%for
disp(sprintf('\ndone!'));

% Output
varargout{1} = atombook;
return;

% --------------------
% subroutines
% --------------------
function np = showprog(m,n,op)
% inputs:
%   m -  current position
%   n -  the total range
%   op - old progress
% outputs:
%   np - new progress

% set parameters
s = 10;         % the step in percent

thres = floor(op+s);
a = floor(m/n*100);
if a >= thres,
    fprintf('* ');
    np = thres;
else
    np = op;
end%if

return;
