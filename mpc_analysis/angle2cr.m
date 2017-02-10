function cr=angle2cr(angle,t,f)
% ANGLE2CR -- convert an angle in the t-f plane to a chirp rate
%
%  Syntax:
%    cr = angle2cr(angle,t,f)
%
%  Inputs:
%    angle  angle (degrees)
%    t      time range
%    f      frequency range
%
%  Outputs
%    cr     chirp rate

% Author(s):    R.J. Cui
% Address:  Sensory communication lab
%           IBBME, U of Toronto
%           4 Taddle Creek Road
%           Toronto, ON M5S 3G9
%           Canada
% Email:    richard.cui@utoronto.ca

% (C)2005 Richard J. Cui; 02-Mar-2005
% $Revision: 0.1 $  $Date: 02-Mar-2005 22:25:16$

narginchk(3,3);

cr = f/t*tan(angle*2*pi/360);
