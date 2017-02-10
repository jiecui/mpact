function varargout = chirpmpr(varargin)

% CHIRPMPR Matching pursuit reconstruction using Gaussian chirplet atoms
%
% [SIG,RTFD] = CHIRPMPR(BOOK,TFDFS,DRAWTFD)
%
%
% Inputs:
%   1. BOOK     - the BOOK structure
%   2. TFDFS    - sampling frequency for time-frequency distribution reconstruction
%   3. DRAWTFD  - to draw the reconstructed time-frequency distribution
%
% Outputs:
%   1. SIG      - reconstructed signal
%   2. RTFD     - rconstructed Wigner-Ville distribution of the signal

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C) 2004-2005 Richard J. Cui, 1-1-2004, created
% $Revision: 0.7 $  $Date: 09-Apr-2005 21:38:52$

% Parsing inputs
book    = varargin{1};
tfdfs   = varargin{2};
drawtfd = varargin{3};

% Parameter setting
L = book.size;
N = book.N;         % N must be an even number
i0 = book.i0;
a = book.rad;
rtfd = zeros(N*tfdfs/2,N*tfdfs);
n = (0:N-1)';
% w = (0:N/2-1)';   % 'w' has not been use yet.
tfn = (0:tfdfs*N-1)'/tfdfs;
tfw = (0:tfdfs*N/2-1)'/tfdfs;

% Error checking

% Main body
sig = zeros(N,1);

for l = 1:L,
    [Rf0gbetal gbetal] = readbook(book,l);
    k = gbetal(1); m = gbetal(2);
    g = gkmn(k,m,i0,a,N);
    % time and frequency shift
    q = gbetal(3); p = gbetal(4);
    % reconstruct the signal
    gq = circshift(g,q);
    phase = exp(sqrt(-1)*2*pi/N*p*n);
    gqp = gq.*phase;
    sig = sig+Rf0gbetal*gqp;
    % Reconstruct the W-V distribution
    alpham = getalpham(k,m,i0,a);
    s = a^k;
    nn = repmat(tfn,1,N*tfdfs/2)';
    ww = repmat(tfw,1,N*tfdfs);
    T = (((nn-q)*cos(alpham))+(ww-p)*sin(alpham))/s;
    W = s*(-(nn-q)*sin(alpham)+(ww-p)*cos(alpham));
    wgbetal = wg(T,W,N);
%     rtfd = rtfd+wgbetal;
    rtfd = rtfd+abs(Rf0gbetal)*wgbetal;        % To weight the signal.
end%for

% Draw TFD if required
if drawtfd
    xid = [1,N]; yid = [1,N/2];
    figure; set(gcf,'NumberTitle','off', 'Name','Wigner-Ville Distribution');
    imagesc(xid,yid,rtfd);colormap(1-gray);axis xy;
    title('Time-frequency visualization of constructive code');
    xlabel('Time');ylabel('Frequency');
    drawnow;
end%if

% Outputs
varargout{1} = sig;
varargout{2} = rtfd;

return;

%-------------------------------------------------------------------
function varargout = wg(varargin)

% TFD = WG(T,W,N)
% Sampling rate normalized to 1/sqrt(2*pi/N)

T = varargin{1};
W = varargin{2};
N = varargin{3};

tfd = sqrt(2/pi)*exp(-2*pi/N*(T.^2+W.^2));

varargout{1} = tfd;

return;