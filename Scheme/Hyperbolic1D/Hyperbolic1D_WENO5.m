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

% Lax-Friedrichs split
FU2 = 0.5*(a*U + alpha*U); % positive
FU1 = 0.5*(a*U - alpha*U); % negative

% Reconstruction of f(u)
[q2,q1] = WENO5_reconstruction(FU2,FU1);

% right-hand side of  du/dt = L(u)
df2 = (q2-circshift(q2,[0,1]))/dx;
df1 = (q1-circshift(q1,[0,1]))/dx;
Lu = -(df2 + df1);