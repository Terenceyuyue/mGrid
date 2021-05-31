function Lu = Euler1D_WENO5(U,dx,alpha,gamma)
% Derive the right-hand side of the semi-discrete scheme: du/dt = L(u)
%
% D.E.	dU/dt + dF(U)/dx = 0,  where
%       U = [rho, rho*u, E]'
%       F(U) = [rho*u, rho*u^2+p, u*(E+p)]'
%
% Spatial discretization: 5th order WENO scheme
% Conserved form:
%   dU_j            F_{j+1/2} - F_{j-1/2}
%  ------   =  -  ------------------------ ,  j = 1,2,...,N
%    dt                      hi
%                    hi
%  where F_{j+1/2} is numerical flux.
%
% The Lax-Friedrichs split is used for stability.
%
% Copyright (C)  Terence Yu.

% Lax-Friedrichs split
FU2 = 0.5*(F(U,gamma) + abs(alpha).*U); % positive
FU1 = 0.5*(F(U,gamma) - abs(alpha).*U); % negative

% Reconstruction of F(U)
[q2,q1] = WENO5_reconstruction(FU2,FU1); % FU2 \approx q2 

% Right-hand side of du/dt = L(u)
df2 = (q2-circshift(q2,[0 +1]))/dx; % circshift: x_{i-1/2}
df1 = (q1-circshift(q1,[0 +1]))/dx; % circshift: x_{i-1/2}
Lu = -(df2 + df1);

end

% Compute flux vector
function flux = F(U,gamma)    
    % F(U) = [rho*u, rho*u^2+p, (E+p)*u] 
    
    % primitive variables
    rho = U(1,:); 
    u   = U(2,:)./rho; 
    E   = U(3,:); 
    p   = (gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector
    flux = [rho.*u; rho.*u.^2+p; (E+p).*u];
end