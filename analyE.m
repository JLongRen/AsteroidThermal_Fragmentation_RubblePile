function T_av=analyE(R,t,t0,XX,phi,rho0)
Cp=837;
lambda=9.5 *10^-7 ;
K=4.3; 
rho=(1-phi)*rho0;
Kappa=K/(rho*Cp);
yr2s=365.25*24*3600;
k=Kappa;
Ho=1.9*10^-7;
To=250;

A1=exp(-lambda*t0);
A=A1*Ho*rho;
fun=cell(XX,1);

%%
p=(lambda/yr2s/k)^0.5;


ts=t-t0;

fun1= @(x) 4*pi.*x.^2*k*A/(K*lambda/yr2s)*exp(-lambda*ts).*((R.*sin(x*p))./(x*sin(R*p))-1);

T1= integral(fun1,0,R)/(4/3*pi*R^3);

for i=1:XX
    
    fun{i,1}= @(x) 2*R^3*A./(x*pi^3*K)*4*pi.*x.^2*(-1)^i/(i*(i^2-lambda/yr2s*R^2/(k*pi^2))).*sin(i*pi*x/R)*exp(-k*i^2*pi^2*ts*yr2s/R^2);
end
for i=1:XX
q(i) = integral(fun{i,1},0,R);
end
T2=sum(q,2)/(4/3*pi*R^3);

T_av=To+T1+T2;

