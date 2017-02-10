function [k,m] = seq2idx(seq, i0, a)
% SEQ2IDX convert sequence of an atom to scale and rotation indexes
% 
% Syntax:
%   [k,m] = seq2idx(seq,i0,a);
%
% Input(s):
%   seq     - sequence number
%   i0      - the first level to rotate
%   a       - the radix of scale
%
% Output(s):
%
% Example:
%
% See also mp_chirplet_bultan.

% Copyright 2004-2016 Richard J. Cui. Created: Mon 01/05/2016 10:33:22.814 AM
% $ Revision: 0.5 $  $ Date: Tue 12/06/2016  6:01:18.628 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

if seq <= i0
    k = seq-1;
    m = 0;
else
    % T(r) = 4*a^(2*r)-1
    % Get k
    d = seq-i0;
    r = 0; findk = 0;
    t1 = 0;
    while ~findk
        t2 = t1+4*a^(2*r)-1;
        if d>t1 && d <=t2
            k = r+i0;
            findk = 1;
        else
            r = r+1;
            t1 = t2;
        end%if
    end%while
    % Get m
    r = 0:k-i0-1;
    t = i0+sum(4*a.^(2*r)-1);
    m = seq-t-2*a^(2*(k-i0));
end % if

end % funciton

% [EOF]