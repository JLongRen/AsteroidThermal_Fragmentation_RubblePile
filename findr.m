function [T,c]=findr(xfL,xfR,R_0,T_0)

x1=xfL-R_0;
x1(x1<0) = nan;
[r1,a]=min(x1);

x2=R_0-xfR;
x2(x2<0) = nan;
[r2,b]=min(x2);
c=b-a-1; %c is how many T point inside the grid
if c==0
    T=T_0(a);
else
    rshell(1)=xfL;
     rshell(2:c+1)=R_0(a+1:a+c);
     rshell(c+2)=xfR;
     T_shell(:)= T_0(a:a+c);
    for i=1:c
      V=4/3*pi*(rshell(2:c+2).^3-rshell(1:c+1).^3);
       T=V*T_shell'/sum(V);
    end
end

