function wvd = chirplet_spgrm(P, n, f, varargin)
% CHIRPLET_SPGRM construct chirplet spectrogram
%
% Syntax:
%   wvd = chirplet_spgrm(P, n, f)
%   wvd = chirplet_spgrm(P, n, f, 'Method', method_name)
% 
% Input(s):
%   P       - matrix of parameters [amp time freq chirp_rate duration; ...]
%   n       - vector of time range (assume fs = 1 Hz)
%   f       - vector of frequency range ([0 1] range)
%   Method  - name-value pair where method_name = {'explicit' 'direct'}.
% 
% Output(s):
%   wvd     - WVD of chirplet spectrogram
% 
% Example:
%
% Note:
%
% References:
%
% See also .

% Copyright 2016 Richard J. Cui. Created: Tue 12/13/2016  2:44:02.637 PM
% $Revision: 0.1 $  $Date: Tue 12/13/2016  2:44:02.638 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
validmethod = {'explicit' 'direct'};

p = inputParser;
p.addRequired('P', @isnumeric);
p.addRequired('n', @isnumeric);
p.addRequired('f', @isnumeric);
p.addParameter('Method', 'explicit', ...
    @(x)any(validatestring(lower(x), validmethod)));

p.parse(P, n, f, varargin{:});
inputs = p.Results;

% =========================================================================
% choose method and estimate spectrogram
% =========================================================================
switch inputs.Method
    case 'explicit'
        wvd = chrpltwvd_explicit(P, n, f);
    case 'direct'
        wvd = chrpltwvd_direct(P, n, f);
end % switch

end % function chirplet_spgrm

% [EOF]
