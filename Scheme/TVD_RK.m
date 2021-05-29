function u = TVD_RK(Lfun,u0,dt)

u1 = u0 + dt*Lfun(u0);
u2 = 3/4*u0 + 1/4*u1 + 1/4*dt*Lfun(u1);
u = 1/3*u0 + 2/3*u2 + 2/3*dt*Lfun(u2);