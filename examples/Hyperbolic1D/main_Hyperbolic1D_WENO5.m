%% main_Hyperbolic1D_WENO5.m (FDM)
% Note: It can also be viewed as a FVM for 1-D problems
%
% Solve the initial boundary value problem of the first-order hyperbolic equation
%
%  D.E.   u_t + a*u_x = 0  ( a>0 ),  x \in (xL, xR), t>0
%  I.C.   u(x,0) = u0(x)
%  B.C.   periodic
%
% The exact solution u(x,t) = u0(x) is periodic w.r.t. time.
%
% The D.E. will be rewritten as
%
%       dU/dt + dF(U)/dx = 0,  where
%       U = u,  F(U) = a*u
%
% Mesh: xL = x0 < x1 < ... < x_N = xR
% cell: I_j = (x_{j-1/2}, x_{j+1/2}),  x_j is the midpoint of I_j
%
% Spatial discretization: 5th order WENO
% Stencil: S0 = { I_{j-2}, I_{j-1}, I_j  }
%          S1 = { I_{j-1}, I_j,     I_{j+1}  }
%          S2 = { I_{j,    I_{j+1}, I_{j+2}  }
% Note: the outside stencils are given by periodic extension
% The semi-discrete conserved scheme is
%
%   du_i            f_{i+1/2} - f_{j-1/2}
%  ------   =  -  ------------------------ ,  i = 1,2,...,N
%    dt                      hi
%
%  which can be written as an ODE system
%
%     du/dt = L(u),   where
%     u = [u1,u2,...,uN]'
%     L(u) is given by WENO discretization
%
% Temporal discretization:
% third order TVD Runge-Kutta method for du/dt = L(u)
%
%    u1 = u^n + dt*L(u^n)
%    u2 = 3/4*u^n + 1/4*u1 + 1/4*dt*L(u1)
%    u^{n+1} = 1/3*u^n + 2/3*u2 + 2/3*dt*L(u2)
%
% References:
%  [1] https://github.com/wme7/WENO
%  [2] Jiang, Guang-Shan, Shu, Chi-Wang. Efficient implementation of weighted
%      ENO schemes. J. Comput. Phys. 126 (1996), no. 1, 202-228.
%
% Copyright (C)  Terence Yu

clc; clear; close all;

%% Parameters
a = 1;
t0 = 0;  tf = 2.0;  
xL = -1; xR = 1;

Nx = 400; % number of cells
dx = (xR-xL)/Nx;
x = xL+dx/2:dx:xR;   % i = 1:N (delete 0)

CFL = 0.4;
dt = CFL*dx/abs(a);

%% initial value
IC = 4;
u_init = @(x) hyperbolic1D_data(x,IC);
u0 = u_init(x);

%% Finite difference method
% right-hand side of the semi-discrete scheme: du/dt = L(u)
Lfun = @(u) Hyperbolic1D_WENO5(u,a,dx);
it = 0;
for t = dt:dt:tf    
    % third order TVD Runge-Kutta method
    u = TVD_RK(Lfun,u0,dt);
    % update
    u0 = u;
    it = it+1;
    % show
    if rem(it,10) == 0
        plot(x,u_init(x),'-k',x,u,'--*','linewidth',1);
        xlim([xL,xR]);
        title('Exact solution is periodic w.r.t. time');
        drawnow;
    end
end