%% main_iMHD1d_WENO5.m (FDM)
% Note: It can also be viewed as a FVM for 1-D problems
%
% The equations will be rewritten as
%
%       dU/dt + dF(U)/dx = 0
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

clc; clear; close all;

%% PDE data
% test cases
IC = 1; % 1-13
pde = iMHD1d_data(IC);
Nx = pde.Nx; dx = pde.dx; x = pde.x;
tf = pde.tf;
cfl = pde.cfl;
gamma = pde.gamma;
vx = pde.vx;
cf = pde.cf;
U0 = pde.U;

% time step
alpha = max(abs(vx)+cf); % using the largest eigenvalue
dt = cfl*dx/alpha;

% boundary condition of periodic extension
% is determined by exact solution
bc = [1:3, Nx-2:Nx];  % stencil width = 3

%% Finite difference method
t = 0;
Bx = pde.Bx;
while t<tf
    % current time
    t = t+dt;
    if t>tf; t = tf; end

    % right-hand side of du/dt = L(u)
    Lfun = @(u) iMHD1d_WENO5(u,dx,alpha,gamma,Bx);

    % third order TVD Runge-Kutta method
    U1 = U0 + dt*Lfun(U0);
    U2 = 3/4*U0 + 1/4*U1 + 1/4*dt*Lfun(U1);
    U = 1/3*U0 + 2/3*U2 + 2/3*dt*Lfun(U2);
    U(:,bc) = pde.U(:,bc);

    % conserved variables
    rho = U(1,:);
    Mx = U(2,:);
    My = U(3,:);
    Mz = U(4,:);
    E = U(5,:);
    By = U(6,:);
    Bz = U(7,:);
    % primitive variables
    vx = Mx./rho;
    vy = My./rho;
    vz = Mz./rho;
    v2 = vx.^2+vy.^2+vz.^2;
    B2 = Bx.^2+By.^2+Bz.^2;
    p = (gamma-1)*(E-0.5*rho.*v2-0.5*B2);
    a2 = gamma*p./rho;	% square of sound speed
    b2 = B2./rho;
    ca2 = Bx.^2./rho;  % square of Alfven speed
    cf = sqrt(0.5*(a2+b2+sqrt((a2+b2).^2-4*a2.*ca2))); % fast speed

    % update
    alpha = max(abs(vx)+cf);
    dt = cfl*dx/alpha;
    U0 = U;
end

figure,
subplot(3,3,1),
plot(x,rho,'--sr'); ylabel('\rho');
subplot(3,3,2),
plot(x,p,'--sr'); ylabel('P');
subplot(3,3,3),
plot(x,E,'--sr'); ylabel('E');
subplot(3,3,4),
plot(x,vx,'--sr');  ylabel('v_x');
subplot(3,3,5),
plot(x,vy,'--sr');  ylabel('v_y');
subplot(3,3,6),
plot(x,vz,'--sr');  ylabel('v_z');
subplot(3,3,7),
plot(x,By,'--sr');  ylabel('By');
subplot(3,3,8),
plot(x,Bz,'--sr');  ylabel('Bz');
subplot(3,3,9),
By(abs(By)<1e-4) = 1; % 0/0 --> 0/1
psi = atand(Bz./By); % Inverse tangent, result in degrees
plot(x,psi,'--sr');  ylabel('\psi');
ylim([min(psi)-0.1, max(psi)+0.1])