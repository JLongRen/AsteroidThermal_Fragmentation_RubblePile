
clear;
clc;
close all;

%% General Physical Properties

R=100*10^3; %unit: m
V=4/3*pi*R^3;
Nt=1500;
H0=1.9 *10^-7 ; %W kg^-3
Cp=837;         %J kg^-1 K^-1
lambda=9.5 *10^-7 ;  %yr^-1
K=4.3;      % W m^-1 K^-1
rho=3720;  % H chondrites (Consolmagno et al. 2008)
M=V*rho;    %unit: kg 
yr2s=365.25*24*3600;  %s/yr

% initial and ending conditions
To=250; % unit:K ambient T
tf=2.27*10^6; %formation time CAI yr
te=tf+250e6; %yr
 
Nx1=100;
Nx2=100;
To_0=To*ones(Nx1,1);

Tas=273.15;  %K
%% Fragmentation Time
 tor=union(linspace(tf,te,Nt),linspace(tf,tf+1e6,100));
 Lt0=size(tor,2);
    
 [ur,grcx0,grdx0,~]=coolFd(R,To_0,tor,Nx1,zeros(Nx1,1),rho);  
 tEs=zeros(1,Nt);

 for i=1:Nt
   tEs(i)=total_energy(grcx0,ur(:,i),rho,Cp,grdx0);
 end

 Tsmean=tEs/(M*Cp);   %average temperature K
 [Tmax1,tmax]=max(ur(1,:)-Tas);

 tBn=tmax; %break up at the peak central T time step 
 tz=tor(tmax); %fragmentation time in CAI

%%  Fagments mass distribution 
det=1;  % unit: year, reassembly time interval
t10=yr2s*det;
Nt1=100;
t1=union(linspace(0,det,Nt1),linspace(0,1,2*Nt1).^3*det);   %linspace(0,det,Nt1) yrs
Lt=size(t1,2);
b=0.99;
Z=200;  %free to change 

phi_f=0.1;   % porosity of all fragments
rho2=rho*(1-phi_f);
Kappa2=K*exp(-phi_f/0.08)/(rho2*Cp);

Repsilon=@(dett, epsilon) sqrt(Kappa2*pi^2*dett*yr2s/log(6/(pi^2*epsilon)));
req=Repsilon(det,0.001);           %m
em=(1-phi_f)*rho*4/3*pi*req.^3;   %kg

[Rs,Ms,m1,mk,Rp,Mp,Np,Y,massleft]=Nff(Z,M,b,rho,0.001,0.9,em);
mk1=1-mk;
Mk=M*massleft;
Mk2=M*mk1;

%% Fragments radial position

R_T=zeros(1,Z+Y+2);
R_T(1)=(Ms(1)*3/(4*pi*rho))^(1/3);
Ms0=Ms(1);
for i=2:Z
 Ms0=Ms0+Ms(i);
 R_T(i)=(Ms0*3/(4*pi*rho))^(1/3);   
end
for i=1:Y
Ms0=Ms0+Mp(i)*Np(i);
R_T(Z+i)=(Ms0*3/(4*pi*rho))^(1/3);   
end
Ms0=Ms0+Mk;
Ms0/M   %check mass conservation

R_T(Z+Y+1)=(Ms0*3/(4*pi*rho))^(1/3);   
R_T(2:Z+Y+2)=R_T(1:Z+Y+1);
R_T(1)=0;

T_R=analyEr(R_T,tz,tf,R,0,20,rho);  %absolute temperature from analytical solution

%% Fragments initial T

Rs1=Rs*(1/(1-phi_f))^(1/3);
Rp1=Rp*(1/(1-phi_f))^(1/3);

tE1=zeros(1,Z);
tEY=zeros(1,Y);

dT=zeros(Nx1,Lt,Z);
d2T=zeros(Nx2,Lt,Y);
u2=zeros(Nx1,Lt,Z);
usmall=zeros(Nx2,Lt,Y);
Tc_r=T_R-To;
xD=linspace(0.005,0.995,100);


%% Fragmentation Cooling 

% Total energy for a cooling solid sphere with uniform initial T
   usum2=@(t)  1-6*sqrt(t)/sqrt(pi)+3*t;  
   for i=1:20    %20 can be changed to 30, 50, 100,...
        usum2= @(t) usum2(t)-12*sqrt(t).*(1/sqrt(pi)*exp(-i^2./t)-i./sqrt(t).*erfc(i./sqrt(t)));
   end

% first 200 largest fragment
for i=1:Z  
    T3=analyef(T_R(i),Rs1(i),det*yr2s,20,0.1,rho);
    tE1(i)=T3*Ms(i)*Cp;
end

% the rest unequilibrium fragments
for i=1:Y

    tc=Rp1(i)^2/Kappa2;
    u=usum2(t10/tc);
    u1=u*Tc_r(Z+i)+To;

    tEY(i)=u1*Mp(i)*Cp;   %total energy for binned fragments
end

tE0=To*(massleft*M)*Cp;    %fragments small enough with ambient T
tE2=sum(tE1)+sum(Np*tEY')+tE0; %add energy from To to total energy(Z+Np*Mp*Cp)

% arrange T in the order of fragments (not in a regualr radial coordinate)
T_rea(1:Z)=tE1./(Cp*Ms);
T_rea(Z+1:Z+Y)=tEY./(Cp*Mp);
T_rea(Z+Y+1)=To;


%%  assign T to the reassembled rubble pile

Nxz=400;
nn1=200;
phi_i=0.1;
phi_t=0.2;
phi(1:nn1)=phi_i; phi(nn1+1:Nxz)=phi_t;
aT3=zeros(1,Nxz);
R_p=2*(3*V/(4*pi*((1-phi_i)+7*(1-phi_t))))^(1/3);


Grid.xmin = 0; Grid.xmax =R_p; Grid.Nx =Nxz; %
Grid = build_grid(Grid); 
grcx2=Grid.xc;

xfL=zeros(1,Nxz);
xfR=zeros(1,Nxz);
cc=zeros(1,Nxz);
for i=1:nn1
    xfL(i)=(1-phi_i)^(1/3)*Grid.xf(i);
    xfR(i)=(1-phi_i)^(1/3)*Grid.xf(i+1);
   [aT3(i),cc(i)]=findr(xfL(i),xfR(i),R_T,T_rea);
end

for i=nn1+1:Nxz
      xfL(i)=((1-phi_i)*(0.5*R_p)^3+ (1-phi_t)*(Grid.xf(i)^3-(0.5*R_p)^3))^(1/3);
      xfR(i)=((1-phi_i)*(0.5*R_p)^3+ (1-phi_t)*(Grid.xf(i+1)^3-(0.5*R_p)^3))^(1/3);
  [aT3(i),cc(i)]=findr(xfL(i),xfR(i),R_T,T_rea);
end
aT3(Nxz)=To;


%%  comparsion before and after fragmentation

figure(3)
hold on
plot(grcx0/1000,ur(:,tmax)-Tas,'r','LineWidth',1)
plot(grcx2/1000,aT3-Tas,'b','LineWidth',1)
plot(grcx0/1000*R_p/R,ur(:,tmax)-Tas,'y','LineWidth',1)
hold off

grid on
legend('before fragmentation','after fragmentation','streched initial body')
xlabel('R (km)')
ylabel('T (^{\circ}C)')

%% Energy Calculation

DtE2=sum(tE2)-To*M*Cp;
T3=analyE(R,tz,tf,20,0,rho); %total energy increase since tf for the parent body
tE3=(T3-To)*M*Cp;

RR=DtE2/tE3     % energy ratio: after/before fragmentation

%% rubble pile

aTz=aT3;
Ntk=500;
dtzz=1e5;
tos=linspace(tz+det,te,Nt);
toz=union(linspace(tz+det,tz+det+dtzz,Ntk),tos);
Ltz=size(toz,2);
% rubble pile profile: u5
[u5,~,~,Krc]=coolFd(R_p,aTz',toz,Nxz,phi',rho);   


%% Discretization for fragment dimensionless solution in the dT-T space

% set up radius and cooling time of the fragment
t_diment=20*yr2s;
Rz=1e3;  % unit: m chooes what radius of fragment you want here 
tc=Rz.^2/Kappa2;
tsD=t_diment./tc;

Nt1=100;
Rsa=1;  
Tc=850+Tas-To;

det2=1e-4;
t1=union(linspace(0,det,Nt1),linspace(0,0.19,3*Nt1).^3*det2);  
t1=union(t1,linspace(0,1,Nt1).^3*0.001*det);
t1=union(t1,linspace(0,1,Nt1).^3*0.1*det);
t3=union(logspace(-7,-1,20),linspace(1e-2,1e-1,10));
t3=union(t3,tsD);
t1=union(t1,t3);

nt3=size(t3,2);
[m3,ntyr] = find(t1==tsD);
    Nx1=1000;
    rs=Rsa;
    rz=1e-4;
    rmp=0.01*linspace(1,99,200);
    Grid.xmin = 0; Grid.xmax =rs; Grid.Nx =Nx1; %
    Grid = build_grid(Grid); 
    
    rx0=union(linspace(0,rz,100),rmp);
    rx=rs-rx0;
    r1=union(Grid.xc,rx);
        r1=union(r1,linspace(0.99,0.99995,100)*rs);
        rsample=1-[1e-4 1e-3 1e-2 0.1 0.5 0.9];
        r1=union(r1,rsample);  
        rsampleT=[0.997 0.99 0.96 0.9 0.67 0.5];
        r1=union(r1,rsampleT);
        
        
   %%  dimensionless solution

   [XD,TD]=meshgrid(r1,t1);
   [u,du]=coolAd(XD,TD,1,rs,20);
   u60=u';
    dT60=du';   %temperature change per 0.1 tc
    rmz=size(r1,2);    

    u6=u60;
    dT6=dT60;
    u6(dT60>-1e-6) = nan;
    dT6(dT60>-1e-6) = nan;
    u6(:,1) = nan;

    N2 =6; % divide dT-T into N2 areas with different transparents
    nr_sample=zeros(1,N2);
    nrnx=zeros(1,N2);

    for i=1:N2
        [mr,nr_sample(i)] = find(r1==rsample(i));
    end

    for i=1:N2
        [mr2,nrnx(i)] = find(r1==rsampleT(i));
    end
    nrnx=[nrnx 1];
 
%% thermal data

Myr=10^6;
%cooling rates======================================
%REE
s1=[839	813	838	824	829	864]; %REE
s2=[10	10	10	10	10	0.5];
%TCa-in-Ol
T_cal=[728	717	712	706	712	733	];
dT_cal=[0.23 0.051 0.079 0.034	0.14 0.12];
% 244Pu Fission Track 
Pux=[200 200 200];
Puy=[2.55 2.95 2.65]/Myr;
ynP = [0.35 0.45 0.35]/Myr; %error data are not shown
ypP = [0.35 0.45 0.35]/Myr;
% Pb-Pb to Ar-Ar 380
PAIx=[380 380 380 380];
PAIy=[7.41 9.80 8.51 7.41]/Myr;
%Metallographic (500 deg C)
sLowT=500*ones(8,1);
sLowdT=[10 15 15 10 10 50 25 30]/Myr;
% Ol-Spinel (700-850 deg C) 
OSG1=775;  
OSG2=75/Myr;
%Opx-Cpx Profile
OCpx=[850 807];
OCpy=[0.1 0.2];

%Chronological Data======================================
%Ar-Ar
Ar_Ary=[280, 280, 280, 240, 240, 240, 240, 240];
Ar_Arx=[103.5 84.5 69.5 58.5 108.5 88.5 118.5 108.5];
yneg = [20 20 20 120 120 120 120 120 ];
ypos = [20 20 20 120 120 120 120 120 ];
xneg = [5 6 6 2.3 1.2 1.5 1.1 1.2 ];
xpos = [5 6 6 2.3 1.2 1.5 1.1 1.2 ];
%Pb-Pb
Pby=[480 450 450 540];
Pbx=[77 64.5 46.5 62.6];
yneg2 = [50 50 50 0];
ypos2 = [50 50 50 0];
xneg2 = [16 1 1 1.4 ];
xpos2 = [16 1 1 1.4 ];
%Hf-W
Hf_Wy=[875 875 837.5 837.5];
Hf_Wx=[10 9.4 12.5 9.5];
yneg3 = [75 75 37.5 37.5];
ypos3 = [75 75 37.5 37.5];
xneg3 = [1.7 1.1 1.2 0.8];
xpos3 = [1.7 1.1 1.2 0.8];


%%  final plot
% color and line setting

[col] = marc_colors();
linew=1;
nrn=nr_sample;

figW = 1200; %scw/2;
gapL = 80;  %gap on left boundary
gapR = 40;  %gap on right boundary
gapT = 30;  %gap on top boundary
gapB = 50;  %gap on bottom boundary
gapH = 90;
gapW=90;

subW = (figW-gapL-gapR-gapW)/2;
subH = 0.7*subW;
figH = subH+gapT+gapB;

subx_a =  gapL/figW;
suby_a = gapB/figH;

subx_b =  (gapL+subW+gapW)/figW;
suby_b = suby_a;

f6=figure(6);

set(f6,'Position', [0, 0, figW, figH])

h1=subplot('position',[subx_a suby_a subW/figW subH/figH]);

mk=220;
ersz=5;  elinw=0.8;
plot(tor(1:tBn)/Myr,ur(1,1:tBn)-Tas,'color',col.green,'LineWidth',1)
hold on 
plot(toz/Myr,u5(1,:)-Tas,'color','k','LineWidth',1)
plot(toz/Myr,u5(mk,:)-Tas,'--','color','k','LineWidth',1)
plot(toz/Myr,u5(373,:)-Tas,'-.','color','k','LineWidth',1)
errorbar(Ar_Arx,Ar_Ary,yneg,ypos,xneg,xpos,'s','MarkerSize',ersz,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',elinw,'color',col.green)
errorbar(Pbx,Pby,yneg2,ypos2,xneg2,xpos2,'s','MarkerSize',ersz,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',elinw,'color',col.blue)
errorbar(Hf_Wx,Hf_Wy,yneg3,ypos3,xneg3,xpos3,'s','MarkerSize',ersz,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',elinw,'color',col.purple)
hold off
xlabel('$t$ (Myr)','interpreter','latex')
ylabel('$T \ \mathrm{(^{\circ}C)}$','interpreter','latex')
axis([0 260 -50 1000])
text(0,1050*1.15/1.1-50,'(a)','FontSize',12);
legend({'center of the parent body','center of the rubble pile','$0.53R_{\mathrm{r}}$','$0.9R_{\mathrm{r}}$',...
    'Ar-Ar','Pb-Pb','Hf-W'},'Location','northeast','NumColumns',2,'FontSize',12,'interpreter','latex')
ax = gca;
ax.FontSize = 12; 


h2=subplot('position',[subx_b suby_b subW/figW subH/figH]);

etaT=@(T) erfinv((T-To)/Tc);  %eto_clo in ren2022
dT_pln=@(T,t) etaT(T).* Tc/(t*sqrt(pi)).*exp(-etaT(T).^2);  
Txx=linspace(To,To+Tc-1,200);


u6z=To+u6*Tc-Tas;
dT6z=yr2s*dT6*Tc;

Snrnx=size(nrnx,2);
nzn=5;

for j=1:Snrnx
    A1=u6z(nrnx(j),:);
    C1=-dT6z(nrnx(j),:);
    BB{j} = A1(~isnan(A1));
    DD{j} = C1(~isnan(A1));
end


ia=linspace(Snrnx,2,Snrnx-1); ib=ia-1;

for j=1:Snrnx-1
    XX{j}=[BB{ia(j)} flip(BB{ib(j)})];
    YY{j}=[DD{ia(j)} flip(DD{ib(j)})];

end

Snrnx2=Snrnx;
transD=linspace(1,0,Snrnx2);
transDm=0.5*(transD(1:end-1)+transD(2:end));
p0=semilogy(s1(1),s2(1),'ks','LineStyle','none');%REE

hold on

dTdc=yr2s*Tc./tc;
semilogy(To+u6(:,ntyr)*Tc-Tas,dTdc*(-dT6(:,ntyr)),'r','LineWidth',1.5);  %from analytical solution

for i=1:Snrnx2-1
pfill=fill(XX{i},YY{i}/tc,col.blue,'FaceAlpha',transDm(i));
pfill.LineStyle = 'none';
end

for i=1:Snrnx2
ALPn{i,1}=[num2str(r1(nrnx(i)),3)];
end

colorticks=flip(transD);
colorLabels=flip(ALPn');
colorLabels{1,1}=0;
colormap(gray(256));
c=colorbar('Ticks',colorticks,'TickLabels',colorLabels);
c.Label.String = '$r''$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
c.FontSize = 12; 
c.FontSize = 12; 

dT_pln1=dT_pln(Txx,yr2s);
dT_pln2=dT_pln(Txx,20*yr2s);
dT_pln3=dT_pln(Txx,1000*yr2s);

semilogy(Txx-Tas,yr2s*dT_pln1,'k','Linewidth',0.8);  % from linear approximation
semilogy(Txx-Tas,yr2s*dT_pln2,'k','Linewidth',1);
semilogy(Txx-Tas,yr2s*dT_pln3,'k','Linewidth',1.5);

p1=semilogy(s1,s2,'ks','MarkerFaceColor',col.orange);%REE
p2=semilogy(T_cal,dT_cal,'kd','MarkerFaceColor',col.orange); %TCa-in-Ol
p4=semilogy(OSG1,OSG2,'d','color',col.orange,'LineWidth',1); % Ol-Spinel
p3=semilogy(OCpx,OCpy,'+','color',col.orange,'LineWidth',1.5);%Opx-Cpx Profile
p6=semilogy(PAIx,PAIy,'o','color',col.orange); %Pb-Pb to Ar-Ar
p5=semilogy(sLowT,sLowdT,'*','color',col.orange); %metallographic
p7=semilogy(Pux,Puy,'^','color',col.orange); %24Pu Fission

TLmid=400;
text(TLmid,6e2,'1 yr')
text(TLmid,3e1,'20 yr')
text(TLmid,5e-1,'1000 yr')

text(700,3*max(dTdc*(-dT6(:,ntyr))),[num2str(t_diment/yr2s,3) ' yr'],'Color','r')

hold off
axis([-100 900 1e-6 1e7])
xlabel('$T \ \mathrm{(^{\circ}C)}$','interpreter','latex')
ylabel('$\dot{T} \ \mathrm{(^{\circ}C/yr)}$','interpreter','latex')

text(-100,LdeL1(-6,7,1.15/1.1),'(b)','FontSize',12);

legend([p1 p2 p3 p4 p5 p6 p7],{'REE','Ca-Ol','Opx-Cpx',...
  'Ol-Spinel','Metallographic','Pb-Pb to Ar-Ar','$^{244}\mathrm{Pu}$ Fission track'},'NumColumns',2,'FontSize',11,'interpreter','latex');
ax = gca;
ax.FontSize = 12; 


