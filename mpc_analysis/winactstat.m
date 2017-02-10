function snrseg = winactstat(x,P,L)
% WINACTSTAT -- Statistics of the results of winACT
%  Syntax:
%    SNRSEG = WINACTSTAT(X,P,L)
%
%  Inputs:
%    X  -   original signal
%    P  -   the estiamted chirplet parameters of each segments
%           p = [m x 5] where m is the number of total segments
%           for each segment, the parameters are [(complex)amplitude,
%           time-center (in samples), frequency-center (in rad),
%           chirp rate (in rad/sample), time-spread (samples)]
%    L  -   the size of each segment
% 
%   Outputs:
%    snrseg -   SNR of each segments
%
% Reconstruct the signal from the results of winact. Note the parameters 
% must be transformed accroding to the chirplet definition in the thesis.

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; Created 03-Nov-2005
% $Revision: 0.1 $  $Date: 03-Nov-2005 09:13:27$

% Parameter setting
M = size(P,1);          % the number of segments
snrseg = zeros(M,1);    % SNR
% To make the parameter consistent with Jeff's notation
Pj = P;
Pj(:,4) = Pj(:,4)*2;
Pj(:,5) = Pj(:,5)/sqrt(2);
% Pj(:,5) = Pj(:,5)*sqrt(2);

% Analysis
for k = 1:M,
    sk = real(make_chirplets(L,Pj(k,:)));
    xk = x((k-1)*L+1:k*L);
    rk = xk - sk;
    snrk = 10*log10(var(sk)/var(rk));
    snrseg(k) = snrk;
end%for

end%function