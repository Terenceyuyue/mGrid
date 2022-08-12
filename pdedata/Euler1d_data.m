function pde = Euler1d_data(IC,x,t)
% x: row vector,  t: scalar
%

switch IC
    case 1
        pde = sinxt(x,t);
        tf = 2;   cfl = 0.5;        
    case 2  % Sod's problem
        pde = Euler_init(x,t,'Sod'); 
        tf = 0.25; cfl = 0.9;        
    case 3  % Lax problem
        a = -0.5;  b = 0.5;
        x = (b-a)/(x(end)-x(1))*(x-a) + a;
        pde = Euler_init(x,t,'Lax'); 
        tf = 0.16; cfl = 0.5;
    case 4  % Right expansion and left strong shock
        pde = Euler_init(x,t,'Right-Left');
        tf = 0.1; cfl = 0.9;
    case 5  
        pde = Euler_init(x,t,'Double-expansion');
        tf = 0.03; cfl = 0.9;
    case 6  
        pde = Euler_init(x,t,'Mach3');
        tf = 0.09; cfl = 0.9;
    case 7  
        pde = Euler_init(x,t,'Shocktube-supersonic');
        tf = 0.162; cfl = 0.9;
    case 8  
        pde = Euler_init(x,t,'Shu-Osher');
        tf = 0.1; cfl = 0.5;
end
pde.tf = tf;  pde.cfl = cfl;

end

%% Case 1: rho = 1 + 0.2*sin(x-t)
function pde = sinxt(x,t)
% x: row vector, t: scalar

% constants
gamma = 1.4;

% primitive variables
rho = 1 + 0.2*sin(x-t);
u = 1 + 0*x;
p = 1 + 0*x;
e =  p./((gamma-1).*rho); % internal energy
E = p./((gamma-1))+0.5*rho.*u.^2;  % total energy density
s =  sqrt(gamma*p./rho);	% speed of sound

% conserved variable
U = [rho; rho.*u; E];

% struct
pde = struct('rho',rho, 'u',u, 'p',p, 'e',e, 'E',e, 's',s, 'U',U);
pde.gamma = gamma;

end  

%% Case 2: Sod
function pde = Euler_init(x,t,ProbName)

% constant
gamma = 1.4; % adiabatic index

% ----------------------- initial values -----------------------
% primitive variables
rho0 = zeros(size(x)); 
u0 = zeros(size(x)); 
p0 = zeros(size(x));

xm = (x(1)+x(end))/2;
id1 = x<=xm;   id2 = ~id1;
switch ProbName
    case 'Sod'        
        rho0(id1) = 1;       rho0(id2) = 0.125;
        u0(id1)   = 0;       u0(id2) = 0;
        p0(id1)   = 1;       p0(id2) = 0.1;
    case 'Lax'
        rho0(id1) = 0.445;    rho0(id2) = 0.5;
        u0(id1)   = 0.698;    u0(id2)   = 0;
        p0(id1)   = 3.528;    p0(id2)   = 0.571; 
    case 'Right-Left'
        rho0(id1) = 1;        rho0(id2) = 1;
        u0(id1)   = 0;        u0(id2)   = 0;
        p0(id1)   = 7;        p0(id2)   = 10;  
    case 'Double-expansion'
        rho0(id1) = 1;        rho0(id2) = 2.5;
        u0(id1)   = -2;       u0(id2)   = 2;
        p0(id1)   = 40;       p0(id2)   = 40;  
    case 'Mach3'
        rho0(id1) = 3.857;    rho0(id2) = 1;
        u0(id1)   = 0.92;     u0(id2)   = 3.55;
        p0(id1)   = 10.333;   p0(id2)   = 1;  
    case 'Shocktube-supersonic'
        rho0(id1) = 1;        rho0(id2) = 0.02;
        u0(id1)   = 0;        u0(id2)   = 0;
        p0(id1)   = 1;        p0(id2)   = 0.02;  
    case 'Shu-Osher'
        rho0(id1) = 3.85714;  rho0(id2) = 1 + 0.2*sin(50*x(id2)-25);
        u0(id1)   = 2.629369; u0(id2)   = 0;
        p0(id1)   = 10.333333;p0(id2)   = 1;         
end


% conserved or other related variables
e0 = p0./((gamma-1).*rho0); % internal energy
E0 = p0./((gamma-1))+0.5*rho0.*u0.^2;  % total energy density
U0 = [rho0; rho0.*u0; E0];  % vector of conserved quantities
s0 = sqrt(gamma*p0./rho0);	% speed of sound

% struct
if t==0
    pde = struct('rho',rho0, 'u',u0, 'p',p0, 'e',e0, 'E',E0, 's',s0, 'U',U0);
    pde.gamma = gamma;
    return;
end

% t>0
% ---------------- exact solution at time t -----------------
[rho,u,p,e,E,s,U] = EulerExact(x,rho0,u0,p0,t,gamma);
pde = struct('rho',rho, 'u',u, 'p',p, 'e',e, 'E',E, 's',s, 'U',U);
pde.gamma = gamma;

end

%% Exact Riemann Solver
function [rho,u,p,e,E,s,U] = EulerExact(x,rho0,u0,p0,tf,gamma)
% Classical Gas Exact Riemann Solver for solving shock-tube problems
% Coded by Manuel Diaz, IAM, NTU 03.09.2011.
%
% This programs was modified by Manuel Diaz, and is based on the code of
% [1] P. Wesseling. PRINCIPLES OF COMPUTATIONAL FLUID DYNAMICS
%     Springer-Verlag, Berlin etc., 2001. ISBN 3-540-67853-0
%     See http://dutita0.twi.tudelft.nl/nw/users/wesseling/
%
% NOTE:
%     A Cavitation Check is the is incorporated in the code. It further
%     prevents plotting for possible but physically unlikely case of
%     expansion shocks.

% Problem definition: Conditions at time t=0
%   rho1, u1, p1
%   rho4, u4, p4
% 'tf' and 'n' are the final solution time and the gas DoFs.

rho1 = rho0(1);   rho4 = rho0(end);
u1 = u0(1);       u4 = u0(end);
p1 = p0(1);       p4 = p0(end);

alpha = (gamma+1)/(gamma-1);

% Assumed structure of exact solution
%
%    \         /      |con |       |s|
%     \   f   /       |tact|       |h|
% left \  a  /  state |disc| state |o| right
% state \ n /    2    |cont|   3   |c| state
%   1    \ /          |tinu|       |k|   4
%         |           |ity |       | |

PRL = p4/p1;
cright = sqrt(gamma*p4/rho4);
cleft  = sqrt(gamma*p1/rho1);
CRL = cright/cleft;
MACHLEFT = (u1-u4)/cleft;

% Basic shock tube relation equation (10.51)
f = @(P) (1+MACHLEFT*(gamma-1)/2-(gamma-1)*CRL*(P-1)/sqrt(2*gamma*(gamma-1 ...
    +(gamma+1)*P)))^(2*gamma/(gamma-1))/P-PRL;

% solve for P = p34 = p3/p4
p34 = fzero(f,3);

p3 = p34*p4;
rho3 = rho4*(1+alpha*p34)/(alpha+p34);
rho2 = rho1*(p34*p4/p1)^(1/gamma);
u2 = u1-u4+(2/(gamma-1))*cleft*(1-(p34*p4/p1)^((gamma-1)/(2*gamma)));
c2 = sqrt(gamma*p3/rho2);
spos = 0.5+tf*cright*sqrt((gamma-1)/(2*gamma)+(gamma+1)/(2*gamma)*p34)+tf*u4;

x0 = 0.5;
conpos = x0 + u2*tf+tf*u4;	% Position of contact discontinuity
pos1 = x0 + (u1-cleft)*tf;	% Start of expansion fan
pos2 = x0 + (u2+u4-c2)*tf;	% End of expansion fan

% primitive variables
p = zeros(size(x));
u= zeros(size(x));
rho = zeros(size(x));
Mach = zeros(size(x));
cexact = zeros(size(x));

for i = 1:length(x)
    if x(i) <= pos1
        p(i) = p1;
        rho(i) = rho1;
        u(i) = u1;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = u(i)/cexact(i);
    elseif x(i) <= pos2
        p(i) = p1*(1+(pos1-x(i))/(cleft*alpha*tf))^(2*gamma/(gamma-1));
        rho(i) = rho1*(1+(pos1-x(i))/(cleft*alpha*tf))^(2/(gamma-1));
        u(i) = u1 + (2/(gamma+1))*(x(i)-pos1)/tf;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = u(i)/cexact(i);
    elseif x(i) <= conpos
        p(i) = p3;
        rho(i) = rho2;
        u(i) = u2+u4;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = u(i)/cexact(i);
    elseif x(i) <= spos
        p(i) = p3;
        rho(i) = rho3;
        u(i) = u2+u4;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = u(i)/cexact(i);
    else
        p(i) = p4;
        rho(i) = rho4;
        u(i) = u4;
        cexact(i) = sqrt(gamma*p(i)/rho(i));
        Mach(i) = u(i)/cexact(i);
    end
end

% conserved and other related quantities
e = p./((gamma-1).*rho); % internal energy
E = p./((gamma-1))+0.5*rho.*u.^2;  % total energy density
U = [rho; rho.*u; E];  % vector of conserved quantities
s = sqrt(gamma*p./rho);	% speed of sound
end