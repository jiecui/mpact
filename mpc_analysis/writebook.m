function varargout = writebook(varargin)
% WRITEBOOK Write the code into the position of the code book
%
% NEWBOOK = WRITEBOOK(OLDBOOK,RFOCODE,CODE,CORRATIO)
%
% Inputs:
% 1. BOOK    - old code book
% 2. RFOCODE - coplex coeffecient of the code
% 3. CODE    - new code added
% 4. CORRATIO - correlation ratio
% 
% Outputs:
% 1. BOOK    - the updated code book  
%
% Update: Feb. 11, 2004
% Update Jan. 6, 2004
% Richard J. Cui
% Email: richard.cui@utoronto.ca
% 

% Parsing inputs
book = varargin{1};
rf0code = varargin{2};
code  = varargin{3};
cr = varargin{4};

% Error checking

% Main body
curp = book.curpoint;
book.code(curp,:) = [rf0code,code];
book.corratio(curp) = cr;
book.curpoint = book.curpoint+1;        % increase pointer

% Outputs
varargout{1} = book;

end % funciton

% [EOF]