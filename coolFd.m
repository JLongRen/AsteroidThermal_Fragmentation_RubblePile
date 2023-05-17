function [uc,r1,d1x,Krc]=coolFd(R,To,t1,Nx,phi,rho)
% R: radius cot:time interval of cooling ts:absolute time 
%To: initail tmeperature

theta=0.5;
Grid.xmin = 0; Grid.xmax =R; Grid.Nx =Nx; %
Grid = build_grid(Grid); 
[D,G,I]=build_ops(Grid); 
 Cp=837;


 H0=1.9 *10^-7 ;
 lambda=9.5 *10^-7 ;
 Ta=250; %ambient T
 yr2s=365.25*24*3600;
 Lt=size(t1,2);
 Param.dof_dir = Grid.dof_xmax; % identify cells on Dirichlet bnd 
 Param.dof_f_dir =Grid.dof_f_xmax; % identify faces on Dirichlet bnd 
 Param.g = Ta; 
 [B,N,fn] = build_bnd(Param,Grid,I); 

 M=spdiags((1-phi)*rho,0,Nx,Nx);
 un=To;
 uc=zeros(Nx,Lt);
uc(:,1)=un;
     %Kappa=7.6029e-07.*rho.*(1-phi).*(1-phi.^0.333);
     %Kappa=(-30*phi+4.4)/Cp;
     Kappa=4.3*exp(-phi/0.08)/Cp; %this Kappa=Conductivity/Capacity
     Kd=comp_mean0(Kappa,-1,Grid);
     L = -D*Kd*G;
 for i=1:Lt-1
    % Kappa=(Aa+Ba./un).*rho.*(1-phi).*(1-phi.^0.333);
     dt=t1(i+1)-t1(i);
 [Lim,Lex] = build_t_ops(Grid,M,L,dt*yr2s);
 IM=Lim(theta);
 EX=Lex(theta);
  Fs0=EX*un; 
  fs=rho*(1-phi)*H0/(Cp*lambda)*(exp(-lambda*t1(i))-exp(-lambda*t1(i+1)));
  u = solve_lbvp(IM,Fs0+fn+yr2s*fs,B,Param.g,N);  
  un=u;
  uc(:,i+1)=u;
  Krc(:,i)=Kappa;
 end
 r1=Grid.xc;
 d1x=Grid.dx;