function Lu = iMHD1d_WENO5(U,dx,alpha,gamma,Bx)

FU = Flux_iMHD1d(U,gamma,Bx);
Lu = LaxFriedrichs(U,FU,dx,alpha);

end

function flux = Flux_iMHD1d(U,gamma,Bx)    
    
    % conserved variables
    rho = U(1,:); 
    Mx = U(2,:);  
    My = U(3,:); 
    Mz = U(4,:);
    E = U(5,:); 
    By = U(6,:);
    Bz = U(7,:);
    % primitive variables
    vx = Mx./rho;
    vy = My./rho;
    vz = Mz./rho;
    v2 = vx.^2+vy.^2+vz.^2;
    B2 = Bx.^2+By.^2+Bz.^2;
    p = (gamma-1)*(E-0.5*rho.*v2-0.5*B2);
    ps = p + 0.5*B2;
    Bv = Bx.*vx + By.*vy + Bz.*vz;
    
    % flux vector
    flux = [rho.*vx;  
        rho.*vx.^2 + p + B2./2 - Bx.^2; 
        rho.*vx.*vy-Bx.*By; 
        rho.*vx.*vz-Bx.*Bz; 
        (E+ps).*vx-Bv.*Bx; 
        By.*vx-Bx.*vy; 
        Bz.*vx-Bx.*vz];
end






