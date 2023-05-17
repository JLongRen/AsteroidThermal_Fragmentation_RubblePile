function T_av=analyef(u,R,t,XX,phi,rho0)
Cp=837;
%K=4.3;
K=4.3*exp(-phi/0.08);
rho=(1-phi)*rho0;
Kappa=K/(rho*Cp);
To=250;

fun=cell(XX,1);
for i=0:XX
fun{i+1,1}= @(x) 4*pi.*x.^2*1./x.*(erfc(((2*i+1)-x)./(2*(Kappa*t).^0.5))-erfc(((2*i+1)+x)./(2*(Kappa*t).^0.5)));
end
%%

for i=0:XX
q(i+1) = integral(fun{i+1,1},0,R);
end
T2=sum(q,2)/(4/3*pi*R^3);
T_av=u-T2;

