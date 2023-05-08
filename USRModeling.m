%% -----☆------☆----- USR Complex Model's Datasets -----☆------☆-----%%
clear
clc
global USR;
global Results;

% Save Datasets Flag 
ResultFlag=0;% Don't save

%% ------☆------☆------☆----- USR Basic Vars -----☆-----☆------☆--------%%

USR.N=10;
USR.L=0.07;
USR.m=1;
USR.J=0.0016;
USR.ct=0.5;
USR.cn=3;
% USR.ct=0.2639;
% USR.cn=8.4;
USR.RiverVx=0.05;
USR.RiverVy=0.05;
lambda1=2.298810^(-7);
lambda2=4.310310^(-4);
lambda3=2.262910^(-5);
mu=0.3958;

%%  ------☆------☆------☆----Calculate Basic Matrix-----☆-----☆------☆-------%%

USR.A=zeros(USR.N-1,USR.N);
    for ii=1:USR.N-1
        USR.A(ii,ii)=1;
        USR.A(ii,ii+1)=1;
    end
USR.D=zeros(USR.N-1,USR.N);
    for ii=1:USR.N-1
        USR.D(ii,ii)=1;
        USR.D(ii,ii+1)=-1;
    end
USR.e=ones(USR.N,1);
USR.e_=ones(USR.N-1,1);
USR.E=[USR.e,zeros(USR.N,1);zeros(USR.N,1),USR.e];
USR.b=zeros(USR.N-1,1);USR.b(USR.N-1)=1;
% USR.b=ones(USR.N-1,1);
USR.H=zeros(USR.N,USR.N-1);
    for ii=1:USR.N-1
        for jj=ii:USR.N-1
            USR.H(ii,jj)=1;
        end
    end
USR.IN=eye(USR.N);
USR.V=USR.A'*(inv(USR.D*USR.D'))*USR.A;
USR.K=USR.A'*(inv(USR.D*USR.D'))*USR.D;

USR.mu=mu*USR.IN;
USR.LAMBDA1=lambda1*USR.IN;
USR.LAMBDA2=lambda2*USR.IN;
USR.LAMBDA3=lambda3*USR.IN;

%%  ------☆------☆------☆----USR Gait Pattern-----☆-----☆------☆-------%%

USR.tStep=.02;
tEnd=20;
t=0:USR.tStep:tEnd;
tNum=tEnd/USR.tStep+1;

%Random values of gait pattern
alphaOrg=15/180*pi+20/180*pi*rand(1,tNum);
delta=50/180*pi;
USR.delta=50/180*pi;
omegaOrg=60/180*pi+5/180*pi*rand(1,tNum);
phi0Org=0+5/180*pi*rand(1,tNum);


alpha=zeros(size(alphaOrg,1),size(alphaOrg,2));
omega=zeros(size(omegaOrg,1),size(omegaOrg,2));
phi0=zeros(size(phi0Org,1),size(phi0Org,2));

% %No Smooth
% alpha=alphaOrg;
% omega=omegaOrg;
% phi0=phi0Org;
% lambda=omega.*t;

% % Fixed Values
% alpha=20/180*pi*ones(1,tNum);
% delta=40/180*pi;
% omega=60/180*pi*ones(1,tNum);
% lambda=omega.*t;
% phi0=0/180*pi*ones(1,tNum);

%%  ------☆------☆------☆---- Datasets Smooth -----☆-----☆------☆-------%%

Tu=0.4;% Greater, more smooth
for ii=1:tNum
    k=floor(t(ii)/Tu)+1;
    alphaDataRef(ii)=alphaOrg(k);
    alpha(ii)=(alphaOrg(k+1)-alphaOrg(k))/Tu*(t(ii)-(k-1)*Tu)+alphaOrg(k);
    omega(ii)=(omegaOrg(k+1)-omegaOrg(k))/Tu*(t(ii)-(k-1)*Tu)+omegaOrg(k);    
    phi0(ii)=(phi0Org(k+1)-phi0Org(k))/Tu*(t(ii)-(k-1)*Tu)+phi0Org(k);  
end
lambda=omega.*t;

figure(1)
hold on
plot(alphaDataRef,'b--','Linewidth',1.5);
plot(alpha,'r','Linewidth',1.5);
legend('平滑设定值','平滑后');
title('数据平滑对比','FontSize',14);
hold off

%% ------☆-------☆-------☆------Initialization -------☆-------☆--------☆-------%%
%derivative of gait pattern
dalpha=zeros(1,tNum);
dalpha(2:tNum)=diff(alpha)/USR.tStep;
ddalpha=zeros(1,tNum);
ddalpha(2:tNum)=diff(dalpha)/USR.tStep;

dlambda=zeros(1,tNum);
dlambda(2:tNum)=diff(lambda)/USR.tStep;
ddlambda=zeros(1,tNum);
ddlambda(2:tNum)=diff(dlambda)/USR.tStep;

dphi0=zeros(1,tNum);
dphi0(2:tNum)=diff(phi0)/USR.tStep;
ddphi0=zeros(1,tNum);
ddphi0(2:tNum)=diff(dphi0)/USR.tStep;

PHI=zeros(USR.N-1,tNum);%每一列是一个时刻的PHI
for ii=1:tNum
    for jj=1:USR.N-1
        PHI(jj,ii)=sin(lambda(ii)+(jj-1)*delta);
    end
end

dPHI=zeros(USR.N-1,tNum);
dPHI(:,2:tNum)=diff(PHI')'/USR.tStep;
ddPHI=zeros(USR.N-1,tNum);
ddPHI(:,2:tNum)=diff(dPHI')'/USR.tStep;

%数据和导数初始化
ALLtheta=zeros(USR.N,tNum);
ALLdtheta=zeros(USR.N,tNum);
ALLddtheta=zeros(USR.N,tNum);

ALLthetaN=ALLtheta(USR.N,:);
ALLdthetaN=zeros(1,tNum);
ALLddthetaN=zeros(1,tNum);

ALLPx=zeros(1,tNum);ALLdPx=zeros(1,tNum);ALLddPx=zeros(1,tNum);
ALLPy=zeros(1,tNum);ALLdPy=zeros(1,tNum);ALLddPy=zeros(1,tNum);

Vt=zeros(1,tNum);
Vn=zeros(1,tNum);
Pt=zeros(1,tNum);
Pn=zeros(1,tNum);

%Set Pcm Initial Value
ALLPx(1)=0;
ALLPy(1)=0;

%Set Gait Pattern Initial Values
% alpha(1)=15/180*pi;
% omega(1)=60/180*pi;
% phi0(1)=0/180*pi;

OutRestrainNum=[];
WarningExtendRestrainNum=0;

%%  ------☆------☆------☆----USR modelling -----☆-----☆------☆-------%%
for m=1:tNum-1
    JointPHI=alpha(m)*PHI(:,m)+phi0(:,m);
    if(max(abs(JointPHI))>=90/180*pi)
        WarningExtendRestrainNum=WarningExtendRestrainNum+1;
        OutRestrainNum=[OutRestrainNum m];
    end
    disp(m);
    USR.Stheta=diag(sin(ALLtheta(:,m)));
    USR.Ctheta=diag(cos(ALLtheta(:,m)));
      
    %Dynamics 
    [USR.Vax,USR.Vay]=Rotation(ALLtheta(:,m),[USR.RiverVx;USR.RiverVy]);%
    [ALLddPx(m+1),ALLddPy(m+1),ALLdPx(m+1),ALLdPy(m+1),ALLPx(m+1),ALLPy(m+1)]=...
    SolvePcm(ALLdtheta(:,m),ALLddtheta(:,m),ALLdPx(m),ALLdPy(m),ALLPx(m),ALLPy(m));                                   
    % Kinetics
    [ALLddthetaN(m+1),ALLdthetaN(m+1),ALLthetaN(m+1),ALLtheta(:,m+1),ALLdtheta(:,m+1),ALLddtheta(:,m+1)]=...
    SolveThetaN(alpha(m),dalpha(m),ddalpha(m),dlambda(m),PHI(:,m),dPHI(:,m),...
    ddPHI(:,m),ddphi0(m),ALLdtheta(:,m),ALLdthetaN(m),ALLthetaN(m),ALLtheta(:,m));

    %Next Time Step
%     ALLdtheta(:,m+1)=(ALLtheta(:,m+1)-ALLtheta(:,m))/USR.tStep;
%     ALLddtheta(:,m+1)=(ALLdtheta(:,m+1)-ALLdtheta(:,m))/USR.tStep;

    Vt(m+1)=cos(ALLthetaN(m+1))*ALLdPx(m+1)+sin(ALLthetaN(m+1))*ALLdPy(m+1);
    Vn(m+1)=-sin(ALLthetaN(m+1))*ALLdPx(m+1)+cos(ALLthetaN(m+1))*ALLdPy(m+1);
    
    Pt(m+1)=Pt(m)+Vt(m+1)*USR.tStep;
    Pn(m+1)=Pn(m)+Vn(m+1)*USR.tStep;
end

%%  -----☆-----☆-----☆---- Save Data (ResultFlag=1) -----☆-----☆------☆-----%%
if ResultFlag
    Results.alpha=alpha(1:end-1)';
    Results.omega=omega(1:end-1)';
    Results.phi0=phi0(1:end-1)';
    Results.ddthetaN=ALLddthetaN';
    Results.Vt=Vt';
    Results.Vn=Vn';
    Results.thetaN=ALLthetaN';
    
    Results.thetaN_Now=Results.thetaN(1:end-1);
    Results.thetaNNextTime=Results.thetaN(2:end);
    Results.ddthetaN_Now=Results.ddthetaN(1:end-1);
    Results.ddthetaNNextTime=Results.ddthetaN(2:end);
    Results.Vt_Now=Results.Vt(1:end-1);
    Results.Vt_NextTime=Results.Vt(2:end);
    Results.Vn_Now=Results.Vn(1:end-1);
    Results.Vn_NextTime=Results.Vn(2:end);

    Results.InputData=[Results.alpha,Results.omega,Results.phi0, Results.thetaN_Now,...
                                  Results.Vt_Now,Results.Vn_Now];
    Results.OutputData=[Results.thetaNNextTime,Results.Vt_NextTime,Results.Vn_NextTime];
    
%     xlswrite( 'Data.xlsx',Results.InputData,1);
%     xlswrite( 'Data.xlsx',Results.OutputData,2);
    
    xlswrite( 'Data1.xlsx',Results.InputData,1);
    xlswrite( 'Data1.xlsx',Results.OutputData,2);
end

%%  -----☆-------☆-------☆------ Plot -------☆-------☆--------☆--------%%
Plotend=tNum;
PlotItv=1;%integer

%Plot vectors
tPlot=1:PlotItv:Plotend;
VtPlot=Vt(1:PlotItv:Plotend);
VnPlot=Vn(1:PlotItv:Plotend);
PxPlot=ALLPx(1,1:1:end);
PyPlot=ALLPy(1,1:1:end);
thetaNPlot=wrapToPi(ALLthetaN);
ddthetaNPlot=ALLddthetaN;
thetaNPlot=thetaNPlot(1:PlotItv:Plotend);
thetaPlot=wrapToPi(ALLtheta(1,:));
thetaPlot=thetaPlot(:,1:PlotItv:Plotend);

%plot(tPlot,ddthetaNPlot,'+','color',[7/8 7/8 7/8]);

%plot
figure(2)
subplot(1,2,1);
plot(tPlot,VtPlot,'r--','Linewidth',1);
legend('\it{V_t(m\cdot s^{-1})}');
xlabel('\it{tStep(s)}');
ylabel('\it{V_t(m\cdot s^{-1})}');
title('\it{USR}切向速度\it{V_t}数据');
subplot(1,2,2);
plot(tPlot,VnPlot,'g--','Linewidth',1);
legend('\it{V_n(m\cdot s^{-1})}');
xlabel('\it{tStep(s)}');
ylabel('\it{V_n(m\cdot s^{-1})}');
title('\it{USR}法向速度\it{V_n}数据');

figure(3)
subplot(2,2,1)
plot(PxPlot,PyPlot,'m','Linewidth',2);
xlabel('\it{P_x(m)}');
ylabel('\it{P_y(m)}');
legend('\it{P_x}和\it{P_y}轨迹');
title('\it{USR}质心位置数据');

subplot(2,2,2);
plot(tPlot,thetaNPlot,'c-','Linewidth',1.5);
legend('\theta _N(rad)');
xlabel('\it{tStep(s)}');
ylabel('\theta _N(rad)');
title('航向角 \it{\theta_N} 数据');

subplot(2,2,3);
plot(tPlot,thetaPlot(1,:),'r-','Linewidth',1.5);
legend('\it{\Theta(rad)}');
xlabel('\it{tStep(s)}');
ylabel('\it{\Theta(rad)}');
title('某节连杆角 \it{\theta} 数据');

figure(4)
    subplot(3,2,1)
    plot(tPlot,alpha,'--','color',[0.2 0.63 0.79]);
    legend('\alpha')
    subplot(3,2,2)
    plot(tPlot,omega,'--','color',[0.2 0.63 0.79]);
    legend('\omega')
    subplot(3,2,3)
    plot(tPlot,phi0,'--','color',[0.2 0.63 0.79]);
    legend('\phi_0')
    subplot(3,2,4)
    plot(tPlot,ALLthetaN,'--','color',[0.2 0.63 0.79]);
    legend('$\bf{\theta_N}$','interpreter','latex')
    subplot(3,2,5)
    plot(tPlot,Vt,'--','color',[0.2 0.63 0.79]);
    legend('V_t')
    subplot(3,2,6)
    plot(tPlot,Vn,'--','color',[0.2 0.63 0.79]);
    legend('V_n')
