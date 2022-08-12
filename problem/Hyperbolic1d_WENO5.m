function Lu = Hyperbolic1d_WENO5(U,a,dx)

% alpha = max |f'(u)|
alpha = abs(a); 

% Lax-Friedrichs splitting
FU = a*U;
Lu = LaxFriedrichs(U,FU,dx,alpha);