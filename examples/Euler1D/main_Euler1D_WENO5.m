%% main_WENO_Euler.m (FDM)
% Note: It can also be viewed as a FVM for 1-D problems
%
% Solve the initial boundary value problem of 1-D Euler system:
%
% D.E.	dU/dt + dF(U)/dx = 0,  where
%       U = [rho, rho*u, E]'
%       F(U) = [rho*u, rho*u^2+p, u*(E+p)]'
%
% I.C.  given by m-file: EulerData.m
% B.C.  periodic
%
% Mesh: xL = x0 < x1 < ... < x_N = xR
% cell: I_j = (x_{j-1/2}, x_{j+1/2}),  x_j is the midpoint of I_j
%
% Spatial discretization: 5th order WENO
%  - Lax-Friedrichs splitting flux is used:
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
%   dU_j            F_{j+1/2} - F_{j-1/2}
%  ------   =  -  ------------------------ ,  j = 1,2,...,N
%    dt                      hj
%
%  which can be written as an ODE system
%
%     dU/dt = L(U),   where
%     U = [U1,U2,...,UN]'
%     L(U) is given by WENO discretization
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
% Copyright (C)  Terence Yu


clc; clear; close all;

%% Mesh
% domain
t0 = 0; tf  = 2;
xL = 0; xR = 1;
% x
Nx = 200; % number of cells
dx = (xR-xL)/Nx;
x = xL+dx:dx:xR;   % j = 1:N (delete 0)
% t
CFL = 0.5;

%% PDE data
% test cases
IC = 2;
pde = Euler1D_data(IC,x,t0,tf);
tf = pde.tf; % use the default value
gamma = pde.gamma;
u0 = pde.u; s0 = pde.s; U0 = pde.U;

% time step
alpha = max(abs(u0)+s0); % using the largest eigenvalue
dt = CFL*dx/alpha;

% boundary condition of periodic extension
% is determined by exact solution
bc = [1:3, Nx-2:Nx];  % stencil width = 3

%% Finite difference method
it = 0; t = 0;
while t<tf
    % current time
    t=t+dt;
    if t>tf; t = tf; end
    
    % right-hand side of du/dt = L(u)
    Lfun = @(u) Euler1D_WENO5(u,dx,alpha,gamma);
    
    % exact solution at time t
    pde = Euler1D_data(IC,x,t,tf);
    
    % third order TVD Runge-Kutta method
    U1 = U0 + dt*Lfun(U0);
    U1(:,bc) = pde.U(:,bc); 
    U2 = 3/4*U0 + 1/4*U1 + 1/4*dt*Lfun(U1);
    U2(:,bc) = pde.U(:,bc); 
    U = 1/3*U0 + 2/3*U2 + 2/3*dt*Lfun(U2);
    U(:,bc) = pde.U(:,bc);    

    % primitive quantities
    rho = U(1,:);
    u = U(2,:)./rho;
    E = U(3,:);
    p = (gamma-1)*(E-0.5*rho.*u.^2);
    s = sqrt(gamma*p./rho);
    
    % update
    alpha = max(abs(u)+s);
    dt = CFL*dx/alpha;
    U0 = U;
    it = it+1;    
    
    % show
    if rem(it,10) == 0
        subplot(1,3,1), 
        plot(x,rho,'--or',x,pde.rho,'k','linewidth',2); 
        subplot(1,3,2), 
        plot(x,u,'--or',x,pde.u,'k','linewidth',2);  
        axis([xL xR min(u)-0.1 max(u)+0.1]);
        subplot(1,3,3), 
        plot(x,p,'--or',x,pde.p,'k','linewidth',2);  
        axis([xL xR min(p)-0.1 max(p)+0.1]);
        drawnow;
    end
end