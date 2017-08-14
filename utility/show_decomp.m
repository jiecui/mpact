function show_decomp(s, P, varargin)
% SHOW_DECOMP Display the results of decomposition on t-f plane
%
% Syntax:
%   show_decomp(s, P)
%   show_decomp(____, fig_name, dbs, nfreq, decf, w)
%   show_decomp(____, 'PType', p_type)
% 
% Input(s):
%   s           - (required) the orignal (real) signal
%   P           - (required) chirplet structure
%   fig_name    - (optional) name of the figures (default fig_name = '')
%   dbs         - (optional) color range of picture in dBs (default = 30)
%   nfreq       - (optional) number of points of frequency range (default =
%                 4 x N, where N is signal length)
%   decf        - (optional) sub-sampling factor in time of the stft, must
%                 be an integer (default = 1)
%   w           - (optional) window width or window, must be an odd number
%                 (default = 5)
%   PType       - (parameter) type of format of P, p_type = {'Cohen', 'Oneill'}
%                 (default p_type = 'Oneill')
%
% Output(s):
%
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2016-2017 Richard J. Cui. Created: Sat 12/17/2016 10:37:41.537 AM
% $Revision: 0.8$  $Date: Tue 07/11/2017  3:07:53.287 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% -------------------------------------------------------------------------
% parse inputs
% -------------------------------------------------------------------------
p = check_inputs(s, P, varargin{:});

s   = p.s;
P   = p.P;
fig_name = p.fig_name;
dbs      = p.dbs; % dbs range in dBs (optional, default is 30)
nfreq    = p.nfreq;
decf     = p.decf; % sub-sampling factor in time of the stft (must be integer)
w        = p.w; % window width or window, must be a odd number
p_type   = p.p_type;


% -------------------------------------------------------------------------
% show
% -------------------------------------------------------------------------
x       = hilbert(s);
N       = length(x);

% c% convert to O'Neill's equation if necessary
% ----------------------------------------------
switch p_type
    case 'oneill'
        fs = 1; 
        P_on = P;
    case 'cohen'
        fs = N;
        P_on = PCo2On(N, P);
end % switch

% Show fft spectrogram
% --------------------
[sx, t, f] = real_spec2(x, fs, nfreq, decf, w);   % spectrum without noise
show_tfd(sx, s, dbs, t, f, fs, 'TFNormType', p_type);
set(gcf, 'Name', 'Spectrogram')

% show the clean chirplet spectrogram
% ---------------------------------------
rs = real(make_chirplets(N, P_on)); % reconstructed signal
t_r = linspace(t(1)*fs, t(2)*fs, 2 * N); % time range
f_r = linspace(f(1)/fs, f(2)/fs, nfreq / 2 + 1); % frequency range

% esitmate WVD of chirplets directly
wig = chirplet_wvd(P_on, t_r, f_r);
show_tfd(wig, rs, dbs, t, f, fs, 'TFNormType', p_type); % dB scale
set(gcf, 'Name', ['WVD - ', fig_name])

% use the direct method to calculate chirplet spectrogram (ACS)
wig = chirplet_spgrm(P_on, t_r, f_r, 'Method', 'direct');
show_tfd(wig, rs, dbs, t, f, fs, 'TFNormType', p_type)
set(gcf, 'Name', ['ACS (direct) - ', fig_name])

% use the explicit method to calculate chirplet spectrogram (ACS)
wig = chirplet_spgrm(P_on, t_r, f_r);
show_tfd(wig, rs, dbs, t, f, fs, 'TFNormType', p_type)
set(gcf, 'Name', ['ACS (explicit) - ', fig_name])

end % function show_decomp

% =========================================================================
% subroutines
% =========================================================================
function wig = chirplet_wvd(P, n, f)
% WVD of chirplets P. Assume O'Neill's equation

fs = 1;
N0 = n(1);
N = n(end) + 1;

d = (f(end) - f(1))/length(f);
nfreq = round(fs / d);
fr = round(f * nfreq / fs + 1); % frequency range

cp = make_chirplets(N, P);
tp = real_wigner1(cp, fs, nfreq);
wig = tp(fr, (2*N0+1):2*N);

end % function

function q = check_inputs(s, P, varargin)

valid_eq = {'oneill', 'cohen'};

p = inputParser;
p.addRequired('s', @isnumeric);
p.addRequired('P', @isnumeric);
p.addOptional('fig_name', '', @ischar);
p.addOptional('dbs', 30, @isnumeric);
p.addOptional('nfreq', 0, @isnumeric);
p.addOptional('decf', 1, @isnumeric);
p.addOptional('w', 5, @isnumberic);
p.addParameter('PType', 'ONeill', @(x) any(validatestring(lower(x), valid_eq)));
p.parse(s, P, varargin{:});

q.s     = p.Results.s;
q.P     = p.Results.P;
q.fig_name  = p.Results.fig_name;
q.dbs   = p.Results.dbs;

if p.Results.nfreq == 0
    N = length(q.s);
    q.nfreq = 4*N;
else
    q.nfreq = p.Results.nfreq;
end % if

q.decf  = p.Results.decf;
q.w     = p.Results.w;
q.p_type    = lower(p.Results.PType);

end % function

% [EOF]
