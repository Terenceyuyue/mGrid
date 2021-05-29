%% main_Hyperbolic1D_LaxWendroff.m (FDM)
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
% Lax-Wendroff scheme (lam = dt/dx):
%
%    u_0^{n+1} = u(0, t_{n+1}) = u0(xL-a*t_{n+1}),                 i = 0,
%    u_i^{n+1} = u_i^n  - a*lam/2*( u_{i+1}^n - u_{i-1}^n ) 
%              - (a*lam)^2/2*( u_{i+1}^n - 2*u_i^n + u_{i-1}^n ),  i = 1:N-1
%
%    at x = xR, i.e., i = N, the missing value u_N^{n+1} is obtained by 
%    using the upwind scheme:
%
%    u_N^{n+1} - u_N^n       u_N^n - u_{N-1}^n
%   ------------------- + a ------------------- = 0,  or
%           dt                     dx
%
%    u_N^{n+1} = u_N^n - a*lam*( u_N^n - u_{N-1}^n )
%
% The CFL condition is : a*lambda <= 1, lambda = dt/dx
%
% Copyright (C)  Terence Yu

clc; clear; close all;

%% Parameters
a = 1; 
t0 = 0; tf = 1; xL = -5; xR = 5;

Nx = 100; 
x = linspace(xL,xR,Nx)'; dx = x(2)-x(1);
lam = 0.5*1/a;  % a*lam = 0.5 < 1
dt = lam*dx;

% initial value
u_init = @(x) (1+0*x).*(x<=0);
%u_init = @(x) (1+0*x).*(x>=0 & x<=1);
%u_init = @(x) sin(pi*x);

%% Finite difference method
u0 = u_init(x);  % t_n 
for t = dt:dt:tf
    u = zeros(size(x));
    % i = 0
    u(1) = u_init(xL-a*t);
    % i = 1:N-1: Lax-Wendroff
    u(2:end-1) = u0(2:end-1) - 0.5*a*lam*( u0(3:end) - u0(1:end-2) ) ...
                + 0.5*(a*lam)^2*( u0(3:end) - 2*u0(2:end-1) + u0(1:end-2) ); 
    % i = N: upwind
    u(end) = u0(end) - a*lam*( u0(end) - u0(end-1) );
    % show
    plot(x,u,'-r', ...
         x,u_init(x-a*t),'--b','linewidth',2);
    xlim([xL, xR]);
    legend('Numerical solution', 'Exact solution');
    drawnow;
    % update
    u0 = u;
end

%% Conclusion
% Lax-Wendroff scheme �ڼ�ϵ㴦�����ԵĲ�ͷ��
% accurate at extrema