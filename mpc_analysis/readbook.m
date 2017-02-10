function varargout = readbook(varargin)

% WRITEBOOK Read the code from a position of the code book
%
% [RF0CODE,CODE] = WRITEBOOK(BOOK,POS)
%
% Inputs:
% 1. BOOK     - the code book
% 2. POS      - the position
% Outputs:
% 1. RF0CODE  - coplex coeffecient of the code
% 2. CODE     - the code
% 
%
% Update Jan. 8, 2004
% Richard J. Cui
% Sensory Communication Lab
% RS422, IBBME
% University of Toronto
% Email: richard.cui@utoronto.ca
% 

% Parsing inputs
book = varargin{1};
pos = varargin{2};

% Error checking

% Main body
book.curpoint = pos;
rf0code = book.code(pos,1);
code = book.code(pos,2:end);

% Outputs
varargout{1} = rf0code;
varargout{2} = code;

return;