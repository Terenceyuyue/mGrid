%% main_Hyperbolic1D_Upwind.m (FDM)
% 
% Solve the initial boundary value problem of the first-order hyperbolic equation 
%  
%     u_t + a*u_x = 0  ( a>0 ),  x \in (xL, xR), t>0
%     u(x,0) = u0(x)         
%     u(0,t) = u0(xL-a*t)    % inflow boundary value condition      
%
% The exact solution is: u(x,t) = u0(x-a*t).
%
% Mesh: xL = x0 < x1 < ... < x_{N-1} < x_N = xR
%
% Spatial discretization: upwind scheme
%
%   du_i(t)        u_i(t) - u_{i-1}(t)
%  --------- + a  --------------------  = 0,   i = 1,2,...,N,
%     dt                   hx
%
% which can be written as an ODE system
%
%     du/dt = L(u),   where    
%     u = [u1,u2,...,uN]'
%     L(u) = -a/h * (u - ud),      ud = [ u_0(t), u1,...,u_{N-1} ]'.
%
% Temporal discretization: 
%
% (1) forward Euler for du/dt = L(u)
%    
%     u^{n+1} = u^n + dt*L(u^n)   
%
% (2) third order TVD Runge-Kutta method for du/dt = L(u)
%
%    u1 = u^n + dt*L(u^n)
%    u2 = 3/4*u^n + 1/4*u1 + 1/4*dt*L(u1)
%    u^{n+1} = 1/3*u^n + 2/3*u2 + 2/3*dt*L(u2)
% 

clc; clear; close all;

%% Parameters
a = 1; 
t0 = 0; tf = 1;  xL = -5; xR = 5;

Nx = 100; 
x = linspace(xL,xR,Nx)';  
hx = x(2)-x(1);
lam = 0.5*1/a;  % CFL: a*lam = 0.5 <= 1
dt = lam*hx;

%% PDE data
% test cases
IC = 1; 
pde = hyperbolic1d_data(IC,a,xL,xR);
uexact = pde.uexact;

%% Spatial discretization
method = 1;  % 1: forward Euler,  2: TVD Runge-Kutta
% du/dt = L(u), u = [u1,u2,...,uN]
Lfun = @(u,t) -a/hx*( u - [uexact(xL,t); u(1:end-1)] );  

%% Temporal discretization 
u0 = uexact(x(2:end),t0);  % delete x0, u0 = [u01,u02,...,u0N]
for t = dt:dt:tf
    switch method
        case 1   % forward Euler
            u = u0 + dt*Lfun(u0,t); 
        case 2   % third order TVD Runge-Kutta method            
            u1 = u0 + dt*Lfun(u0,t);  
            u2 = 3/4*u0 + 1/4*u1 + 1/4*dt*Lfun(u1,t);  
            u = 1/3*u0 + 2/3*u2 + 2/3*dt*Lfun(u2,t);  
    end
    uf = [uexact(xL,t); u]; % add u(0,t)
    % show 
    plot(x,uf,'-r', ...
         x,uexact(x,t),'--b','linewidth',2);
    xlim([xL, xR]);
    legend('Numerical solution', 'Exact solution');
    drawnow; 
    % update
    u0 = u;
end

%% Conclusion
% Upwind scheme is stable
% smooth at discontinuities 
% degenerated accuracy at smooth extrema ( u_0 = @(x) sin(pi*x) )
