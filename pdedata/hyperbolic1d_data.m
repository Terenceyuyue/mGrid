function pde = hyperbolic1d_data(IC,a,xL,xR)
% References:
%  [1] Jiang, Guang-Shan, Shu, Chi-Wang. Efficient implementation of weighted
%      ENO schemes. J. Comput. Phys. 126 (1996), no. 1, 202-228.
%  [2] K. Ou, P. Vincent, and A. Jameson. High-order methods for diffusion
%      equation with energy stable flux reconstruction scheme.
%

uexact = @(x,t) uinit(x-a*t,IC,xL,xR); 
pde.uexact = uexact;

end

function u0 = uinit(x,IC,xL,xR)

Lx = xR-xL; x = x/(0.5*Lx);  % [-a, a] --> [-1 1]

switch IC
    case 1 % Riemann problem
        u0 = (1+0*x).*(x<=0);
    case 2 % Square Jump
        u0 = (1+0*x).*(x<=0.2 & x>=-0.2);
    case 3 % Oleg's trapesoidal
        u0 = exp(0.1-x).*(x<=0.2 & x>=-0.2);
    case 4 % combination of Gaussians ..., see Example 1 in Ref. [1]
        u0 = combination(x);   % [-1 1]
    case 5 % Gaussian wave
        u0 = exp(-20*x.^2); 
    case 6 % sin 
        u0 = sin(pi*x);
    case 7 % hyperbolic tangent
        mu = 0.02;
        u0 = 0.5*(1-tanh(x/(4*mu)));
end

end

function u0 = combination(x)
% The solution contains a smooth but narrow combination of Gaussians,
% a square wave, a sharp triangle wave, and a half ellipse.

% x-ranges
x1 = x.*(x>=-0.8 & x<=-0.6);
x2 = x.*(x>=-0.4 & x<=-0.2);
x3 = x.*(x>= 0.0 & x<= 0.2);
x4 = x.*(x>= 0.4 & x<= 0.6);

% constants
a = 0.5; z = -0.7; delta = 0.005;
alpha = 10; beta = log(2)/(36*delta^2);

% functions
G = @(x,b,r) exp(-b*(x-r).^2);
F = @(x,d,r) sqrt(max(1-d^2*(x-r).^2,0));

% Initial condition
u0 = 1/6*( G(x1,beta,z-delta) + G(x1,beta,z+delta) + 4*G(x1,beta,z)) + ...
    1*(x2~=0) + ...
    1-abs(10*(x3-0.1)) + ...
    1/6*(F(x4,alpha,a-delta) + F(x4,alpha,a+delta) + 4*F(x4,alpha,a));
end
