function Lu = Euler1d_WENO5(U,dx,alpha,gamma)


FU = Flux_Euler1d(U,gamma);
Lu = LaxFriedrichs(U,FU,dx,alpha);

end

function flux = Flux_Euler1d(U,gamma)    
    % F(U) = [rho*u, rho*u^2+p, (E+p)*u] 
    
    % primitive variables
    rho = U(1,:); 
    u   = U(2,:)./rho; 
    E   = U(3,:); 
    p   = (gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector
    flux = [rho.*u; rho.*u.^2+p; (E+p).*u];
end