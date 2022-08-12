%% main_Hyperbolic1D_WENO5.m (FDM)
% Note: It can also be viewed as a FVM for 1-D problems
%
% Solve the initial boundary value problem of the first-order hyperbolic equation
%
%  D.E.   u_t + a*u_x = 0  ( a>0 ),  x \in (xL, xR), t>0
%  I.C.   u(x,0) = u0(x)
%  B.C.   periodic
%
% The exact solution u(x,t) = u0(x-a*t) 
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
%  - Lax-Friedrichs splitting is used:
%
%     F(u) = F+(u) + F-(u)
%     F+(u) = 1/2* ( F(u) + alpha*u )
%     F-(u) = 1/2* ( F(u) - alpha*u )
%     alpha = max |F'(u)|
%
%  - Stencil for positive flux at x_{j+1/2}: 
%          S0+ = { I_{j-2}, I_{j-1}, I_j  }
%          S1+ = { I_{j-1}, I_j,     I_{j+1}  }
%          S2+ = { I_{j,    I_{j+1}, I_{j+2}  }
%  - Stencil for negative flux at x_{j+1/2}: 
%          S0- = { I_{j-1}, I_j,      I_{j+1}  }
%          S1- = { I_j,     I_{j+1},  I_{j+2}  }
%          S2- = { I_{j+1}, I_{j+2},  I_{j+3}  }
%  - The outside stencils are given by periodic extension, in this case the
%    first or the last stencil should be deleted.
%
% The semi-discrete conserved scheme is
%
%   du_j            f_{j+1/2} - f_{j-1/2}
%  ------   =  -  ------------------------ ,  j = 1,2,...,N
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
%  [3] 张德良. 计算流体力学教程. 北京: 高等教育出版社, 2010, pp. 424-430.
%

clc; clear; close all;

%% Parameters
% domain
a = 1;
t0 = 0; tf = 1;  
xL = -5; xR = 5;

% x
Nx = 200; % number of cells
dx = (xR-xL)/Nx;
x = xL+dx:dx:xR;   % j = 1:N (delete 0)
% t
CFL = 0.5;
dt = CFL*dx/abs(a);

%% PDE data
% test cases
IC = 4; % 1-7
pde = hyperbolic1d_data(IC,a,xL,xR);
uexact = pde.uexact;
% boundary condition of periodic extension 
% is determined by exact solution
bc = [1:3, Nx-2:Nx]; 

%% Finite difference method
% right-hand side of the semi-discrete scheme: du/dt = L(u)
u0 = uexact(x,t0);
Lfun = @(u) Hyperbolic1d_WENO5(u,a,dx);
it = 0;
for t = dt:dt:tf    
    % third order TVD Runge-Kutta method
    u1 = u0 + dt*Lfun(u0);
    u1(bc) = uexact(x(bc),t);
    u2 = 3/4*u0 + 1/4*u1 + 1/4*dt*Lfun(u1);
    u2(bc) = uexact(x(bc),t);
    u = 1/3*u0 + 2/3*u2 + 2/3*dt*Lfun(u2);
    u(bc) = uexact(x(bc),t);
    % update
    u0 = u;
    % show
    plot(x,u,'-k', ...
         x,uexact(x,t),'--or','linewidth',2);
    %axis([xL xR min(u)-0.1 max(u)+0.1]);
    legend('Numerical solution', 'Exact solution');
    drawnow;
end