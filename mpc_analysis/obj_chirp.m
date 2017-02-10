function [fval, g] = obj_chirp(P, x)
% OBJ_CHIRP objective function minimized by FINMINCON
% 
% It usually provides with objective gradient and Hessian.
% 
%  Syntax:
%   [fval, g] = obj_chirp(P, x)
%
%  Inputs:
%   P       - estimate of chirplet parameters ([time freq chirp_rate duration])
%   x       - signal being fitted
%
%  Outputs:
%   f       - value of the objective function
%   g       - gradients of f function

% Copyright 2005-2016 Richard J. Cui. Created: Tue 02/22/2005 10:33:22.814 AM
% $ Revision: 0.2 $  $ Date: Tue 12/06/2016  6:01:18.628 PM $
%
% 3236 E Chandler Blvd Unit 2036
% Phoenix, AZ 85048, USA
%
% Email: richard.jie.cui@gmail.com

% FVAL
N = length(x);
fval = -abs(x' * make_chirplets(N, [1 P]))^2;

% Gradient
% f(t,f,c,d) = |z(t,f,c,d)|^2
% df/dt = 2 re{(dz/dt) z*}
t = P(1);       % get the initial estimates
f = P(2);
c = P(3);
d = P(4);

n = (1:N)';
y = make_chirplets(N,[1 P]);
z_conj = conj(sum(x.*conj(y)));

dz_dt = -sum(x.*conj(y).*((n-t)/2/d^2+1i*c*(n-t)+1i*f));
g(1) = 2*real(dz_dt*z_conj);

dz_df = -sum(x.*conj(y).*(-1i*(n-t)));
g(2) = 2*real(dz_df*z_conj);

dz_dc = -sum(x.*conj(y).*(-1i/2*(n-t).^2));
g(3) = 2*real(dz_dc*z_conj);

dz_dd = -sum(x.*conj(y).*((n-t).^2/2/d^3-1/2/d));
g(4) = 2*real(dz_dd*z_conj);

end % function

% [EOF]