function Lu = Hyperbolic1D_WENO5(U,a,dx)
% Derive the right-hand side of the semi-discrete scheme: du/dt = L(u)
%
% D.E. u_t + a*u_x = 0 , a>0 is a constant
%  or  du/dt + dF(u)/dx = 0,  F(u) = a*u
%
% Spatial discretization: 5th order WENO scheme
% Conserved form:
%   du_i            f_{i+1/2} - f_{j-1/2}
%  ------   =  -  ------------------------
%    dt                      hi
%  where f_{i+1/2} is numerical flux.
%
% The Lax-Friedrichs split is used for stability.
%
% Copyright (C)  Terence Yu.

% alpha for Lax-Friedrichs flux
% alpha = max |f'(u)|
alpha = abs(a); 

% positive part
ind = [0 +1]; 
q2 =   WENO5_reconstruction(U,ind);
f2 = 0.5*(a*q2 + alpha*q2);      % at x_{i+1/2}
df2 = (f2-circshift(f2,ind))/dx; % circshift gives the values at x_{i-1/2}

% negative part
ind = [0 -1];
q1 =   WENO5_reconstruction(U,ind);
f1 = 0.5*(a*q1 - alpha*q1);      % at x_{i+1/2}
df1 = (f1-circshift(f1,ind))/dx; % circshift gives the values at x_{i-1/2}

% right-hand side of  du/dt = L(u)
Lu = -(df2 + df1);