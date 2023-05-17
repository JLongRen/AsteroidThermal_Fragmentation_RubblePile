function [u,du]=coolAd(x,t,Kappa,R,N)
u=1;
du=0;

for i=0:N
u=u-R./x.*(erfc(((2*i+1)*R-x)./(2*(Kappa*t).^0.5))-erfc(((2*i+1)*R+x)./(2*(Kappa*t).^0.5)));
end
for i=0:N
    y=((2*i+1)*R-x)./(2*(Kappa*t).^0.5);
    z=((2*i+1)*R+x)./(2*(Kappa*t).^0.5);
du=du-R./(x*sqrt(pi)).*(exp(-y.^2).*y./t-exp(-z.^2).*z./t);
end