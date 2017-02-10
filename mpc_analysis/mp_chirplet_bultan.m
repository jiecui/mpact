function [Rf0gbetal, betal] = mp_chirplet_bultan(x, D, i0, a, verbose)
% MP_CHIRPLET_BULTAN Implements Matching-pursuit with Bultan chirplet atoms
%
% Syntax:
%   [Rf0gbetal, betal] = mp_chirplet_bultan(x, D, i0, a, verbose)
%
% Inputs:
%   x       - the input signal
%   D       - the depth
%   i0      - the first level of rotaiton, typically i0 = 1
%   a       - the radix, typically a = 2
%   verbose - show info
%
% Outputs:
%   Rf0betal- the complex amplitude of the maximum chirplet
%   betal   - the parameter (scale, rotate, time-shift, freq-shift) index of
%             the maximum chirplet
%
% References:
%   [1]	J. Cui, "Adaptive chirplet transform for the analysis of visual 
%       evoked potentials," Doctor of Philosophy Dissertation (Ph.D.),
%       University of Toronto, 2006.
%   [2]	A. Bultan, "A four-parameter atomic decomposition of chirplets,"
%       IEEE Transactions on Signal Processing, vol. 47, pp. 731-745, Mar
%       1999.
%
% Example:
%
% Note:
%   Due to LEM, Mallat & Zhang's fast algorithm won't be employed.
% 
% See also mp_chirplet, seq2idx, writebook.

% Copyright 2016 Richard J. Cui. Created: Wed 12/14/2016  3:59:32.295 PM
% $Revision: 0.1 $  $Date: Wed 12/14/2016  3:59:32.324 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% Parameter setting
% =========================================================================
x = x(:); % force to be column vector
N = length(x); % the length of the signal

r = 0:(D-i0-1); % total # of levels, r: 0,1,...,D-i0-1
nidx = i0+sum(4*a.^(2*r)-1); % total number of atoms, from the level

% =========================================================================
% Now matching the atoms one by one (not the fast algorithm)
% =========================================================================
maxabscoe = 0;              % the max coefficient
perc = 10;                  % percentage finished
% seq: atom sequence (k,m);q: time-shift; p: frequency-shift
if verbose
    fprintf('\n%d chirplet atoms in dictionary.\nSearch finished: ', nidx);
end %if
for seq = 1:nidx           % for each chirplet in the dictionary
    [k,m] = seq2idx(seq,i0,a);          % get index of scale and rotation from the sequence
    g_km = gkmn(k,m,i0,a,N);            % gkmn: Gaussian chirplet atom at scale k and rotation m
    for q = 0:N-1                      % for each time position
        g_kmq = circshift(g_km,q);      % time-shift
        Rf0_gkmq = x.*conj(g_kmq);    % column vector
        Rf0_g = fft(Rf0_gkmq);          % freq by timeshift
        aR = abs(Rf0_g);                % absolute, N by 1
        f_shift = find(aR == max(aR));  % frequency-shift
        p = f_shift(1)-1;               % more than 1, get 1st
        if aR(p+1) > maxabscoe          % find the max
            maxabscoe = aR(p+1);
            Rf0gbetal = Rf0_g(p+1);
            % [k,m] = seq2idx(seq,i0,a);
            betal = [k,m,q,p];          % \beta_l
        end %if
    end %for
    if verbose
        ratio = floor(seq/nidx*100);
        if ratio >= perc
            fprintf(1,'%d%% ',perc);
            perc = perc+10;
        end %if
    end %if
end %for

end % function mp_chirplet_bultan

% [EOF]
