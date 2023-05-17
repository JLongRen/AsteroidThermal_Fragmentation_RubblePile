function   tE=total_energy(rxc,u,rho,Cp,dx)
% author: Jialong Ren 
% date: 12 11 2019
V=4/3*pi*((rxc+0.5*dx).^3-(rxc-0.5*dx).^3); %V(1,N)  u(N,1)
tE=Cp*rho*V*u;

% Grid.xmin = 0; Grid.xmax =2; Grid.Nx =100; %
% Grid = build_grid(Grid); 
% rxc=Grid.xc;
% dx=Grid.dx;
% 
% V=4/3*pi*((rxc+0.5*dx).^3-(rxc-0.5*dx).^3);
% tE=Cp*rho2*V*ones(100,1)*250;
% Tave=tE/(sum(V)*rho2*Cp);