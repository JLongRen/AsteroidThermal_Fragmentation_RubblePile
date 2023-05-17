function [Rs ,Ms,mr,mk,Rp,Mp,Np,Y,massleft,Mdn,Mp0]=Nff(Z,M,b,rho,error,alpha,minm)
mr=1/zeta(1/b);
m_max= M*1/zeta(1/b);
n1=1:Z;
C= m_max^b; 
m=@(n) m_max*n.^-(1/b);
Ms=m(n1);
mk=sum(Ms)/M;
Rs=(3*Ms/(rho*4*pi)).^(1/3);
M0=Ms(Z);
i=1;
massleft(1)=1-mk;
N(1)=Z;
Mp0(1)=M0;
Mz=M0;

Mdn=0;
if massleft<=error
    Mp=0;
Np=0;
Rp=0;
Y=0;
else

while massleft>error && Mz>minm   %minm in kg
    i=i+1;
 Mp0(i)=alpha*Mp0(i-1);  %bottom M

  N(i)=floor(C*Mp0(i)^-b);

Mdn(i)=m_max*(b/(b-1))*(N(i)^(1-1/b)-N(i-1)^(1-1/b));


Mp(i-1)=Mdn(i)/(N(i)-N(i-1));
 Nk(i)=N(i)-N(i-1);
 massleft=massleft-Mdn(i)/M;
 Mz=Mp0(i);
 end

Np(:)=Nk(2:i);
Rp=(3*Mp/(rho*4*pi)).^(1/3);
Y=i-1;
end