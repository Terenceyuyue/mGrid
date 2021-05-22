%% main_Hyperbolic1D_LaxWendroff.m
% 
% Solve the initial value problem of the first-order hyperbolic equation 
%  
%     u_t + a*u_x = 0  ( a>0 )
%     u(x,0) = u0(x)
%
% using the Lax-Wendroff scheme (lam = dt/dx):
%
%    u_i^{n+1} = u_i^n  - a*lam/2*( u_{i+1}^n - u_{i-1}^n ) 
%                       - (a*lam)^2/2*( u_{i+1}^n - 2*u_i^n + u_{i-1}^n )
%
% The CFL condition is : a*lambda <= 1, lambda = dt/dx
% The exact solution is: u(x,t) = u0(x-at)
%
% Note: The missing values of u_i^{n+1} are given by the ones at the last time step
%
% Copyright (C)  Terence Yu

%% Finite difference method
clc; clear; close all;

a = 1; 
t0 = 0; tf = 1;
xL = -5; xR = 5;

Nx = 100; 
x = linspace(xL,xR,Nx)'; dx = x(2)-x(1);
lam = 0.5*1/a;  % a*lam = 0.5 < 1
dt = lam*dx;
u_0 = @(x) (1+0*x).*(x<=0);
%u_0 = @(x) (1+0*x).*(x>=0 & x<=1);

u0 = u_0(x);  % t_n 
i = 1;
for t = t0:dt:tf
    % the missing values are given by u0
    uf = u0;  % t_{n+1}
    % compute by scheme
    uf(i+1:end-1) = u0(i+1:end-1) - 0.5*a*lam*( u0(i+2:end) - u0(i:end-2) ) ...
                  + 0.5*(a*lam)^2*( u0(i+2:end) - 2*u0(i+1:end-1) + u0(i:end-2) );    
    % show
    plot(x,uf,'-r',x,u_0(x-a*t),'--b','linewidth',2);
    legend('Numerical solution', 'Exact solution');
    drawnow;
    % update the location of last one of missing values
    i = i+1;   
    if i>Nx, break; end
    % update
    u0 = uf;
end

%% Conclusion
% Lax-Wendroff 格式稳定
% 在间断点处有明显的波头振荡
