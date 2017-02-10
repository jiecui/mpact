function varargout = getalpham(varargin)
% GETALPHAM Calculate the discrete rotational angles

% angm = getalpham(k,m,i0,a)

% Update: 21-Dec-2005 10:03:34
% 
% Richard J. Cui
% Sensory communication lab
% IBBME, U of Toronto
% 4 Taddle Creek Road
% Toronto, ON M5S 3G9
% Canada
% Email: richard.cui@utoronto.ca
% 

% Parse inputs
k = varargin{1};
m = varargin{2};
i0 = varargin{3};
a = varargin{4};

% Error checking

% Calculation
angm = atan(m/a^(2*(k-i0)));

% Output
varargout{1} = angm;

return;