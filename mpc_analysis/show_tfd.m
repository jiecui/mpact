function show_tfd(tfd, x, varargin)
% SHOW_TFD -- show an image plot of a time-frequency distribution
% 
% Syntax:
%   show_tfd(tfd, x)
%   show_tfd(____,dbs, t, f, fs, fnts, lin)
%   show_tfd(____, 'TFNormType', tf_type)
%
% Inputs:
%   tfd         - time-frequency distribution
%   x           - original signal, assume real signal
%   dbs         - (optional) range in dBs scale (default 25)
%   t           - (optional) range of sampling times (default [1 size(tfd,2)])
%   f           - (optional) range of frequency values (default 0->.5 Hz)
%   fs          - (optional) sampling frequency (default 1 Hz)
%   fnts        - (optional) font size of axis labels (default 10)
%   lin         - (optional) linear scale flag (default false to use dB scale)
%   TFNormType  - (parameter) T-F axis normalization type, 
%                 tf_type = {'', 'Oneill', 'Cohen'}, where 'Oneill' type
%                 normalization is [0 N-1]x[0 \pi], and 'Cohen' type [0
%                 1]x[0 \pi] (defualt '')
% 
% Note:
%   show an image plot of a TFD with a linear or dB amplitude scale
%   together with the original signal and its spectrum for a real signal.
%
% See also .

% Copyright 2005-2017 Richard J. Cui. Created: Tue 03/02/2005 2:46:52.278 PM
% $ Revision: 0.5 $  $ Date: Tue 07/11/2017 12:59:22.887 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% -------------------------------------------------------------------------
% parse inputs
% -------------------------------------------------------------------------
p = check_inputs(tfd, x, varargin{:});
tfd     = p.tfd;
x       = p.x;
dbs     = p.dbs;
t       = p.t;
f       = p.f;
fs      = p.fs;
fnts    = p.fnts;
lin     = p.lin;
tf_type = p.tf_type;

% find the portion of signal x
xa = floor(t(1)*fs)+1;
xb = floor(t(2)*fs)+1;
xseg = x(xa:xb);            % signal segment

% Find the fft of x

nfft = 1024;                % # points of fft
xfft = fft(x,nfft);         % fft of x 
if rem(nfft,2) == 0
    n = nfft/2+1;
else
    n = (nfft+1)/2;
end
famp = abs(xfft(1:n));
fa = floor(f(1)*nfft/fs)+1;
fb = floor(f(2)*nfft/fs)+1;
fampseg = famp(fa:fb);

% plot the image
if lin                              % linear amplitude scale
    mintfd = min(tfd(:));           % normailze between 0 and 1
    tfd = tfd - mintfd;
    maxtfd = max(tfd(:));
    tfd = tfd./maxtfd;
    clim = [0 1];
else                                % dB amplitude scale
    tfd = tfd./sum(sum(tfd));       % normalized the dynamic range to [0 1]
    tfd = 20*log10(abs(tfd)+eps);   % dBs
    tfd = tfd - max(max(tfd));
    clim = [-dbs,0];
end %if

% -------------------------------------------------------------------------
% TF axis normalization
% -------------------------------------------------------------------------
switch tf_type
    case {'cohen', 'oneill'}
        f = f/fs*2;
        x_label = 'Normalized Time (s)';
        y_label = 'Nromalized Frquency (\times\pi rad)';
    otherwise
        x_label = 'Time (s)';
        y_label = 'Frquency (Hz)';        
end % switch

disp_prest(tfd, t, f, xseg, fampseg, clim, fnts, x_label, y_label);

end

% =========================================================================
% subroutines
% =========================================================================
function q = check_inputs(tfd, x, varargin)

valid_eq = {'oneill', 'cohen'};

p = inputParser;
p.addRequired('tfd', @isnumeric);
p.addRequired('x', @isnumeric);
p.addOptional('dbs', 25, @isnumeric);
p.addOptional('t', 0, @isnumeric);
p.addOptional('f', [0, .5], @isnumeric);
p.addOptional('fs', 1, @isnumeric);
p.addOptional('fnts', 10, @isnumeric);
p.addOptional('lin', false, @islogical);
p.addParameter('TFNormType', '', @(x) any(validatestring(lower(x), valid_eq)));
p.parse(tfd, x, varargin{:});

q = p.Results;
q.tf_type = lower(p.Results.TFNormType);
if q.t == 0
    q.t = [1, size(q.tfd, 2)];
end % if

end % function

% [EOF'