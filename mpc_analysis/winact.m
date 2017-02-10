function [P,M] = winact(c,L)
% WINACT -- The windowed ACT method
%  Syntax:
%    [P,M] = WINACT(C,L)
%
%  Inputs:
%    c  -   the signal (complex)
%    L  -   the length of each segment
% 
%  Outputs:
%    P  -   the estiamted chirplet parameters of each segments
%           p = [m x 5] where m is the number of total segments
%           for each segment, the parameters are [(complex)amplitude,
%           time-center (in samples), frequency-center (in rad),
%           chirp rate (in rad/sample), time-spread (samples)]
%           Note: the format is consistent with my thesis
%    M  -   the number of segments
%
% Decompose a signal by using the windowed ACT method.

% Copyright 2005-2016 Richard J. Cui. Created: Tue 02/22/2005 2:46:52.278 PM
% $ Revision: 0.3 $  $ Date: Tue 12/13/2016 10:33:20.464 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% Parameter settings
N = length(c);      % signal size
M = floor(N/L);     % the number of segments
P = zeros(M,5);     % the matrix for the estimated parameters
Q = 1;              % number of atoms to look for
% level = 2;          % level of difficulty (0 ->3)
res = L;            % resolution parameter (optional, default is 64)
verbose = 0;        % verbose flag
% MP parameters
D = 5;          % the depth
i0 = 1;         % first level to rotate
radix = 2;      % radix of scale

% Estmate only ONE chirplet in each segment
for k = 1:M
    ck = c((k-1)*L+1:k*L);  % the data of the kth segment
    Pk = mp_adapt_chirplets(ck,Q,res,D,i0,radix,verbose);  % MP method
    P(k,:) = Pk;
end%for
P = transpara(P);     % trsnsfer to thesis format 

end

% -----------------
% subroutines
% -----------------
function Pout = transpara(Pin)
% to change the estimates according to the equation in the thesis

P = Pin;
P(:,4) = P(:,4)/2;          % change the chirp rate
P(:,5) = P(:,5)*sqrt(2);    % change the time-spread
Pout = P;

end

% [EOF]