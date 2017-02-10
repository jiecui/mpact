function varargout = idx2seq(varargin)

% seq = idx2seq(k,m,i0,a);
%
% look up sequence in the index book
% 1. k - scale index
% 2. m - angle index
% 3. i0 - the first level to rotate
% 4. a - the radix of scale
%
% Update: Jan. 5, 2004
% Richard J. Cui
% Email: richard.cui@utoronto.ca
% 

% Sparse inputs
k = varargin{1};
m = varargin{2};
i0 = varargin{3};
a = varargin{4};

% Error checking

% Calculating
if k<i0
    seq = k+1;
    varargout{1} = seq;
    return;
elseif k ==i0
    t = i0;
else
    r = 0:k-i0-1;
    t = i0+sum(4*a.^(2*r)-1);
end%if
seq = m+t+2*a^(2*(k-i0));

% Ouputs
varargout{1} = seq;

return;