function T_av=analyEr(R_T,t,t0,R,phi,XX,rho0) %%XX is the number of series
Cp=837;
lambda=9.5 *10^-7 ;
K=4.3*exp(-phi/0.08);   %note the K as afunction of phi
rho=(1-phi)*rho0;
Kappa=K/(rho*Cp);
yr2s=365.25*24*3600;
k=Kappa;
Ho=1.9*10^-7;
To=250;
A1=exp(-lambda*t0);
A=A1*Ho*rho;
fun=cell(XX,1);
Lz=size(R_T,2);
%%
p=(lambda/yr2s/k)^0.5;

ts=t-t0;

fun1= @(x) 4*pi.*x.^2*k*A/(K*lambda/yr2s)*exp(-lambda*ts).*((R.*sin(x*p))./(x*sin(R*p))-1);


    for j=1:Lz-1
        T1(j)= integral(fun1,R_T(j),R_T(j+1))/(4/3*pi*(R_T(j+1)^3-R_T(j)^3));
    end
for i=1:XX
      fun{i,1}= @(x) 2*R^3*A./(x*pi^3*K)*4*pi.*x.^2*(-1)^i/(i*(i^2-lambda/yr2s*R^2/(k*pi^2))).*sin(i*pi*x/R)*exp(-k*i^2*pi^2*ts*yr2s/R^2);
end
for i=1:XX
    for j=1:Lz-1
        q(i,j) = integral(fun{i,1},R_T(j),R_T(j+1));
    end
end
   for j=1:Lz-1
        V_R(j) =4/3*pi*(R_T(j+1)^3-R_T(j)^3);
    end
T2=sum(q,1)./V_R;
T_av=To+T1+T2;

