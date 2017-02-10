function [tfd,rsig] = winactrec(P,L)
% WINACTREC -- Reconstruction from the results of WINACT
%  Syntax:
%    [TFD,RSIG] = WINACTREC(P,L)
%
%  Inputs:
%    P  -   the estiamted chirplet parameters of each segments
%           p = [m x 5] where m is the number of total segments
%           for each segment, the parameters are [(complex)amplitude,
%           time-center (in samples), frequency-center (in rad),
%           chirp rate (in rad/sample), time-spread (samples)]
%    L  -   the size of each segment
% 
%   Outputs:
%    tfd    figure of the time-frequency distribution
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

% (C)2005 Richard J. Cui; Created 22-Feb-2005
% $Revision: 0.3 $  $Date: 03-Nov-2005 08:36:42$

% Parameter setting
tfd = [];
M = size(P,1);          % the number of segments
rsig = zeros(1,M*L);    % reconstructed signal
num_t = L*4;            % number of points along time-axis
num_f = L*4;            % number of points along frequency-axis
% To make the parameter consistent with Jeff's notation
Pj = P;
Pj(:,4) = Pj(:,4)*2;
% Pj(:,5) = Pj(:,5)/sqrt(2);
Pj(:,5) = Pj(:,5)*sqrt(2);

% Meshgrid points
[T,W] = meshgrid(linspace(0,L,num_t),linspace(0,pi,num_f));

for k = 1:M,
    % reconsturct signal
    sk = real(make_chirplets(L,Pj(k,:)));
    rsig((k-1)*L+1:k*L) = sk;
    % reconstruct the t-f WVD
    A = P(k,1);         % the complex amplitude
    tc = P(k,2);        % the time-center
    fc = P(k,3);        % the frequency-center
    cr = P(k,4);        % the chirp rate
    ts = P(k,5);        % the time-spread
    wvg = 2*abs(A)*exp(-(T-tc).^2/ts^2)...
        .*exp(-ts^2*((W-fc)-2*cr*(T-tc)).^2);
    tfd = [tfd,wvg];
end%for

end
