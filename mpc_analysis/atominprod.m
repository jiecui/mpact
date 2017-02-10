function varargout = atominprod(varargin)

% ATOMINPROD Calculate inner product of one atom and the dictionary
%
% INPROD = ATOMINPROD(DICT,BETAL,I0,A)
% 
% Inputs:
% 1. DICT
% 2. BETAL
% 3. I0
% 4. A
% 
% Outputs:
% 1. INPROD
%
% Update: Jan. 7, 2004
% Update Jan. 6, 2004
% Richard J. Cui
% Email: richard.cui@utoronto.ca
% 

% Parsing inputs
dict = varargin{1};
betal = varargin{2};
i0 = varargin{3};
a = varargin{4};

% Parameter setting
[N,nidx] = size(dict);
inprod = zeros(nidx,N,N);

% Error Checking

% Main body
kl = betal(1); ml = betal(2); ql = betal(3); pl = betal(4);
l = idx2seq(kl,ml,i0,a);
gklml = dict(:,l);
gklmlql = circshift(gklml,ql);

for q = 0:N-1,
    repgklmlql = repmat(gklmlql,1,nidx);   % repeat gklmlql
    gkmq = conj(circshift(dict,q));
    gg = repgklmlql.*gkmq;
    gginner1 = fft(gg);
    gginner2 = circshift(gginner1,pl); % frequency shift, p x nidx
    inprod(:,q+1,:) = gginner2.';      % must be transpose only, don't do conjugate
end%for

% Outputs
varargout{1} = inprod;
return;