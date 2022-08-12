function pde = iMHD1d_data(IC,Nx)
% References:
%  [1] D. Ryu and T.W. Jones. Numerical magnetohydrodynamics in
%  astrophysics: Algorithm and tests for one-dimensional flow. Astrophys. J. 
%  442: 228-258, 1995.
%  [2] J.M. Stone, T.A. Gardiner, P. Teuben and et al. Athena: A new code
%  for astrophysical MHD. The Astrophysical Journal Supplement Series. 178:
%  137-177, 2008
%

% x: row vector

switch IC
    case 1
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ1a');
        tf = 0.08; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 1a of Ryu & Jones (1995) ----\n')
    case 2
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ1b');
        tf = 0.03; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 1b of Ryu & Jones (1995) ----\n')
    case 3
        if nargin==1, Nx = 1024; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ2a');
        tf = 0.2; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 2a of Ryu & Jones (1995) ----\n')
    case 4
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ2b');
        tf = 0.035; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 2b of Ryu & Jones (1995) ----\n')
    case 5
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ3a');
        tf = 0.01; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 3a of Ryu & Jones (1995) ----\n')
    case 6
        if nargin==1, Nx = 512; end
        xL = -0.5;  xR = 0.5;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ3b');
        tf = 0.1; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 3b of Ryu & Jones (1995) ----\n')
    case 7
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ4a');
        tf = 0.15; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 4a of Ryu & Jones (1995) ----\n')
    case 8
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ4b');
        tf = 0.15; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 4b of Ryu & Jones (1995) ----\n')
    case 9
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ4c');
        tf = 0.15; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 4c of Ryu & Jones (1995) ----\n')
    case 10
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ4d');
        tf = 0.16; cfl = 0.8;
        fprintf('---- problem = Riemann problem from Figure 4d of Ryu & Jones (1995) ----\n')
    case 11 % Brio & Wu shocktube
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ5a');
        tf = 0.1; cfl = 0.8;
        fprintf('---- problem = Brio & Wu shock tube ----\n')
        fprintf('---- /problem = Riemann problem from Figure 5a of Ryu & Jones (1995) ----\n')
    case 12 
        if nargin==1, Nx = 512; end
        xL = 0;  xR = 1;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'RJ5b');
        tf = 0.16; cfl = 0.8;
        fprintf('---- /problem = Riemann problem from Figure 5b of Ryu & Jones (1995) ----\n')
    case 13
        if nargin==1, Nx = 800; end
        xL = -1.5;  xR = 1.5;  dx = (xR-xL)/Nx;
        x = xL+dx:dx:xR;   % j = 1:N (delete 0)
        pde = iMHD_init(x,'Torrilhon'); % Fig. 6 in Torrilhon (2003)
        tf = 0.4; cfl = 0.8;
        fprintf('---- problem = Torrihlon shock tube ----\n')


end
pde.x = x; pde.dx = dx; pde.Nx = Nx;
pde.tf = tf;  pde.cfl = cfl;

end

%% Case 2: Sod
function pde = iMHD_init(x,ProbName)



% ----------------------- initial values -----------------------
% primitive variables
[rho,vx,vy,vz,p,Bx,By,Bz] = deal(zeros(size(x)));

xm = (x(1)+x(end))/2;
id1 = (x<=xm);   id2 = ~id1;
switch ProbName
    case 'RJ1a'
        Bxc = 5/sqrt(4*pi); gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 1;
        vx(id1)  = 10;      vx(id2) = -10;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 20;      p(id2) = 1;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = Bxc;     By(id2) = Bxc;
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ1b'
        Bxc = 3/sqrt(4*pi); gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 0.1;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 1;       p(id2) = 10;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 5/sqrt(4*pi);  By(id2) = 2/sqrt(4*pi);
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ2a'
        Bxc = 2/sqrt(4*pi);          gamma = 5/3;
        rho(id1) = 1.08;             rho(id2) = 1;
        vx(id1)  = 1.2;              vx(id2) = 0;
        vy(id1)  = 0.01;             vy(id2) = 0;
        vz(id1)  = 0.5;              vz(id2) = 0;
        p(id1)   = 0.95;             p(id2) = 1;
        Bx(id1)  = Bxc;              Bx(id2) = Bxc;
        By(id1)  = 3.6/sqrt(4*pi);   By(id2) = 4/sqrt(4*pi);
        Bz(id1)  = 2/sqrt(4*pi);     Bz(id2) = 2/sqrt(4*pi);
    case 'RJ2b'
        Bxc = 3/sqrt(4*pi);  gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 0.1;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 2;
        vz(id1)  = 0;       vz(id2) = 1;
        p(id1)   = 1;       p(id2) = 10;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 6/sqrt(4*pi);  By(id2) = 1/sqrt(4*pi);
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ3a'
        Bxc = 0;            gamma = 5/3;
        rho(id1) = 0.1;     rho(id2) = 0.1;
        vx(id1)  = 50;      vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 0.4;     p(id2) = 0.2;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = -1/sqrt(4*pi);  By(id2) = 1/sqrt(4*pi);
        Bz(id1)  = -2/sqrt(4*pi);  Bz(id2) = 2/sqrt(4*pi);
    case 'RJ3b'
        Bxc = 0;            gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 1;
        vx(id1)  = -1;      vx(id2) = 1;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 1;       p(id2) = 1;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 1;       By(id2) = 1;
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ4a'
        Bxc = 1;            gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 0.2;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 1;       p(id2) = 0.1;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 1;       By(id2) = 0;
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ4b'
        Bxc = 1.3;          gamma = 5/3;
        rho(id1) = 0.4;     rho(id2) = 1;
        vx(id1)  = -0.66991;vx(id2) = 0;
        vy(id1)  = 0.98263; vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 0.52467; p(id2) = 1;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 0.0025293;   By(id2) = 1;
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ4c'
        Bxc = 0.75;         gamma = 5/3;
        rho(id1) = 0.65;    rho(id2) = 1;
        vx(id1)  = 0.667;   vx(id2) = 0.4;
        vy(id1)  = -0.257;  vy(id2) = -0.94;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 0.5;     p(id2) = 0.75;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 0.55;    By(id2) = 0;
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ4d'
        Bxc = 0.7;          gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 0.3;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 1;
        p(id1)   = 1;       p(id2) = 0.2;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 0;       By(id2) = 1;
        Bz(id1)  = 0;       Bz(id2) = 0;
    case 'RJ5a'
        Bxc = 0.75;         gamma = 2;
        rho(id1) = 1;       rho(id2) = 0.125;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 1;       p(id2) = 0.1;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 1;       By(id2) = -1;
        Bz(id1)  = 0;       Bz(id2) = 0;
     case 'RJ5b'
        Bxc = 1.3;         gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 0.4;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 1;       p(id2) = 0.4;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 1;       By(id2) = -1;
        Bz(id1)  = 0;       Bz(id2) = 0;       
    case 'Torrilhon' 
        Bxc = 1;            gamma = 5/3;
        rho(id1) = 1;       rho(id2) = 0.2;
        vx(id1)  = 0;       vx(id2) = 0;
        vy(id1)  = 0;       vy(id2) = 0;
        vz(id1)  = 0;       vz(id2) = 0;
        p(id1)   = 1;       p(id2) = 0.2;
        Bx(id1)  = Bxc;     Bx(id2) = Bxc;
        By(id1)  = 1;       By(id2) = cos(3);
        Bz(id1)  = 0;       Bz(id2) = sin(3);
end

% conserved or other related variables
Mx = rho.*vx;
My = rho.*vy;
Mz = rho.*vz;
v2 = vx.^2+vy.^2+vz.^2;
B2 = Bx.^2+By.^2+Bz.^2;
E = p./(gamma-1) + 0.5*rho.*v2 + 0.5*B2; % total energy density
U = [rho; Mx; My; Mz; E; By; Bz];  % Bx is not included for 1d problems
a2 = gamma*p./rho;	% square of sound speed
b2 = B2./rho;
ca2 = Bx.^2./rho;  % square of Alfven speed
cf = sqrt(0.5*(a2+b2+sqrt((a2+b2).^2-4*a2.*ca2))); % fast speed

% struct
pde = struct('rho',rho, 'vx',vx, 'vy',vy, 'vz',vz, 'p',p, ...
    'Bx',Bx, 'By',By, 'Bz',Bz, 'Mx',Mx, 'My',My, 'Mz',Mz, ...
    'E',E, 'a2',a2, 'b2',b2,'ca2',ca2, 'cf',cf, ...
    'U',U,'v2',v2,'B2',B2);
pde.gamma = gamma;
end % end of iMHD_init