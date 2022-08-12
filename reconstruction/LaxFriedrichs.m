function Lu = LaxFriedrichs(U,FU,dx,alpha)
% Derive the right-hand side of the semi-discrete scheme: du/dt = L(u)
%
% D.E.	dU/dt + dF(U)/dx = 0
% Conserved form:
%
%   dU_j            F_{j+1/2} - F_{j-1/2}
%  ------   =  -  ------------------------ ,  j = 1,2,...,N
%    dt                      hj
%                    
%  where F_{j+1/2} is numerical flux.
%
% The Lax-Friedrichs splitting is used for stability.
%

% Lax-Friedrichs splitting
FU2 = 0.5*(FU + abs(alpha).*U); % positive
FU1 = 0.5*(FU - abs(alpha).*U); % negative

% Reconstruction of F(U)
[q2,q1] = WENO5_reconstruction(FU2,FU1); % WENO5 can be replaced by other methods

% Right-hand side of du/dt = L(u)
df2 = (q2-circshift(q2,[0 +1]))/dx; % circshift: x_{j-1/2}
df1 = (q1-circshift(q1,[0 +1]))/dx; % circshift: x_{j-1/2}
Lu = -(df2 + df1);