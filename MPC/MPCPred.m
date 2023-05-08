clear
clc

%% CpxModel 参数 
global USR;

USR.N=10;
USR.L=0.07;
USR.m=1;
USR.J=0.0016;
USR.ct=0.5;
USR.cn=3;
USR.RiverVx=0.05;
USR.RiverVy=0.05;
lambda1=2.298810^(-7);
lambda2=4.310310^(-4);
lambda3=2.262910^(-5);
mu=0.3958;
USR.delta=50/180*pi;
USR.tStep=.02;

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

%% MPC

import casadi.*
global MPC;
%GPRModelling需要读取全局变量
MPC.ddthetaNModel=load('ddthetaNModel.mat');
MPC.VtModel=load('VtModel.mat');
MPC.VnModel=load('VnModel.mat');

% ddthetaNModel=load('ddthetaNModel.mat');
% VtModel=load('VtModel.mat');
% VnModel=load('VnModel.mat');

%% casADi
%% casADi定义
Np=5;
U_MPC=SX.sym('U',3,Np);
X_MPC=SX.sym('U',3,Np+1);
P_MPC=SX.sym('P',3,1);
X_MPC(:,1)=P_MPC;
N=10;

DefineVarStart=tic;

for ii=1:Np
    X_MPC(:,ii+1)=GPRModelling(X_MPC(:,ii),U_MPC(:,ii),1);
    disp('casADi初始赋值进度：');
    disp(ii);
end

% save('X_MPC.mat','X_MPC');
disp('模型转化为casADi用时:')
disp(toc(DefineVarStart));

%% MPC 构建求解体
SetupNLPStart=tic;

%时间
USR.tStep=0.02;
tEnd=4;
tNum=tEnd/USR.tStep+1;
t=0:USR.tStep:tEnd;
Cnt=1;

%目标函数
ObjFun=0;
gamma=0.025;

for ii=1:Np
    TempU=U_MPC(:,ii);
    TempVt=X_MPC(2,ii);
%     TempVn=abs(X_MPC(1,ii));
%     ObjFun=ObjFun-TempVt+gamma*(TempU')*TempU;%方案1
%     ObjFun=ObjFun-TempVt+TempVn+gamma*(TempU')*TempU;%方案4
%     ObjFun=ObjFun-TempVt;%方案1
%     ObjFun=ObjFun+abs(TempVt-0.04);%方案4
    ObjFun=ObjFun+(TempVt-0.05)^2;%方案4
end

%非线性约束变量 对应下面args
NonLnrCons_PHI=[];
delta=50/180*pi;

for jj=1:N-1
     NonLnrCons_PHI=[NonLnrCons_PHI,U_MPC(1,:).*sin(U_MPC(2,:).^t(Cnt:Cnt+Np-1)+(jj-1)*delta)+U_MPC(3,:)];
end
NonLnrCons_PHI=max(max(abs(NonLnrCons_PHI)));

%NonLnrCons_PHI维度: 1    x   10(N-1)

%NLP求解结构体
OptimU=reshape(U_MPC,3*Np,1);
NlpMPC=struct('f',ObjFun,'x',OptimU,'g',[],'p',P_MPC);

%NLP求解对象
opts = struct;
opts.ipopt.max_iter = 200;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-9;%Acceptable Convergence Tolerance
MPCSolver=nlpsol('MPCSolver', 'ipopt', NlpMPC,opts);

%NLP求解非线性约束
args=struct;
args.lbx=[];
args.ubx=[];
for ii=1:Np
    args.lbx=[args.lbx;15/180*pi;60/180*pi;0/180*pi];
    args.ubx=[args.ubx;35/180*pi;60/180*pi;10/180*pi];
end
% args.lbx=zeros(3*Np,1);
% args.ubx=2*pi*ones(3*Np,1);
% args.lbg=-90/180*pi*ones(1,Np*(N-1));
% args.ubg=90/180*pi*ones(1,Np*(N-1));
args.lbg=[];
args.ubg=[];
disp('构建结构体用时：')
disp(toc(SetupNLPStart));

%% MPC 开始
%定义结果存储变量

XInit=[0;0;0];
% UInit=[15/180*pi;60/180*pi;0/180*pi];
UInit=[0;pi/3;0];
UEveryOptim=zeros(3*Np,1);
UEveryOptim=zeros(3,Np);
UEveryOptim(2,:)=ones(1,Np)*pi/3;
UEveryOptim=reshape(UEveryOptim,[3*Np,1]);

ALLUMPC=zeros(3,tNum-1);
ALLXMPC=zeros(3,tNum);

ALLXMPC(:,1)=XInit;

%Start MPC
StartMPCOptimTimeFlag=tic;

clear 'MPC_CpxModelling';
for m=1:tNum-1
     fprintf('MPC 当前进展：第%d时刻 \n',m);
     args.p=ALLXMPC(:,m);
     args.x0=UEveryOptim;
     OptimAns=MPCSolver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,'lbg', args.lbg,...
                        'ubg', args.ubg,'p',args.p);
     UEveryOptim=full(OptimAns.x);
     ALLUMPC(:,m)=UEveryOptim(1:3);
     ALLXMPC( : ,m+1)=MPC_CpxModelling(ALLXMPC( : ,m),ALLUMPC(:,m));  
%      ALLXMPC( : ,m+1)=GPRModelling(ALLXMPC( : ,m),ALLUMPC(:,m),1);  
      UEveryOptim(1:(end-3))=UEveryOptim(4:end);
      Cnt=Cnt+1;
end
FinishMPCOptimTimeFlag=toc(StartMPCOptimTimeFlag);
fprintf('MPC优化结束，用时%f 秒 \n',FinishMPCOptimTimeFlag);

% for jj=1:N-1
%     PHI(jj,:)=ALLUMPC(1,:).*sin(ALLUMPC(2,:).^t(1:100)+(jj-1)*50/180*pi)+ALLUMPC(3,:);
% end
% MaxPHI=max(max(abs(PHI)));

%% 验证关节角约束
ALLPHI=zeros(N-1,size(ALLUMPC,2));
    for jj=1:N-1
        ALLPHI(jj,:)=ALLUMPC(1,:).*sin(ALLUMPC(2,:).^t(1:end-1)+(jj-1)*delta)+ALLUMPC(3,:);
    end
JudgeMat=(ALLPHI>pi/2)|(ALLPHI<-pi/2);
if(sum(sum(JudgeMat)))
    disp('可惜，最优化后的U序列得到的PHI不满足范围。');
else
    disp('最优化后的U序列得到的PHI均满足范围。');
end

%% MPC 数据结果处理

%存数据标志位
SaveResultFlag=0;
%存入哪个ResultMethod# 文件
SaveFileNum=0;

% ALLXMPC=xlsread('ResultsMethod3.xlsx',2);
% ALLUMPC=xlsread('ResultsMethod3.xlsx',1);

clear 'MPC_CpxModelling';
XANS=ALLXMPC;
for m=1:size(ALLUMPC,2)
   XANS(:,m)=MPC_CpxModelling(ALLXMPC( : ,m),ALLUMPC(:,m));
end
plot(XANS(2,:));

figure(1)
subplot(2,2,1)
plot(ALLXMPC(1,:),'-','color',[0.01,0.66,0.62],'Linewidth',2);
xlabel('TimeStep \it{(0.02s)}','FontSize',14);
legend('$\bf{\ddot{\theta}_N}$','interpreter','latex')
title('$\bf{ MPC: \it{  \theta_N}}$','interpreter','latex','Fontsize',16);
box on
axes('Position',[0.23,0.65 0.2 0.15]);
plot([40:1:100],ALLXMPC(1,40:100),'-','color',[1 0.5 0],'Linewidth',2);
box on

subplot(2,2,2)
plot(ALLXMPC(2,:),'-','color',[0.01,0.66,0.62],'Linewidth',2);
xlabel('TimeStep \it{(0.02s)}','FontSize',14);
legend('V_t','Fontsize',16)
title('\bf{MPC: \it{V_t}}','Fontsize',16);
box on
axes('Position',[0.65,0.63 0.2 0.15]);
plot([40:1:100],ALLXMPC(2,40:100),'-','color',[1 0.5 0],'Linewidth',2);
box on

subplot(2,2,3)
plot(ALLXMPC(3,:),'-','color',[0.01,0.66,0.62],'Linewidth',2);
xlabel('TimeStep \it{(0.02s)}','FontSize',14);
legend('V_n','Fontsize',16)
title('\bf{MPC : \it{V_n}}','Fontsize',16);
box on
axes('Position',[0.23,0.25 0.15 0.15]);
plot([40:1:100],ALLUMPC(3,40:100),'-','color',[1 0.5 0],'Linewidth',2);
box on

figure(2)
subplot(2,2,1)
plot(ALLUMPC(1,2:end),'-','color',[0,0.78,0.55],'Linewidth',1.5);
xlabel('TimeStep \it{(0.02s)}','FontSize',14);
legend(' \alpha','Fontsize',14);
title('\bf{\it{MPC: \alpha}}','FontSize',16);
box on
axes('Position',[0.25,0.77 0.12 0.12]);
plot([30:1:70],ALLUMPC(1,30:70),'-','color',[1 0.5 0],'Linewidth',2);
box on

subplot(2,2,2)
plot(ALLUMPC(2,:),'-','color',[0,0.78,0.55],'Linewidth',1.5);
xlabel('TimeStep \it{(0.02s)}','FontSize',14);
legend(' \omega','Fontsize',14);
title('\bf{\it{MPC: \omega}}','FontSize',16);
box on
axes('Position',[0.65,0.61 0.15 0.12]);
plot([50:1:100],ALLUMPC(2,50:100),'-','color',[1 0.5 0],'Linewidth',2);
box on

subplot(2,2,3)
plot(ALLUMPC(3,3:end),'-','color',[0,0.78,0.55],'Linewidth',1.5);
% axis([0 100 -0.1 0.25]);
xlabel('TimeStep \it{(0.02s)}','FontSize',14);
legend(' \phi_0','Fontsize',14);
title('\bf{\it{MPC: \phi_0}}','FontSize',16);
box on
axes('Position',[0.23,0.25 0.15 0.15]);
plot([50:1:100],ALLUMPC(3,50:100),'-','color',[1 0.5 0],'Linewidth',2);
box on

if SaveResultFlag
    if (SaveFileNum==1)%能耗Vt最小
        xlswrite('ResultsMethod1.xlsx',ALLUMPC,1);
        xlswrite('ResultsMethod1.xlsx',ALLXMPC,2);
    end
    if (SaveFileNum==2)
        xlswrite('ResultsMethod2.xlsx',ALLUMPC,1);
        xlswrite('ResultsMethod2.xlsx',ALLXMPC,2);
    end
    if (SaveFileNum==3)%绝对值趋近于0.08
        xlswrite('ResultsMethod3.xlsx',ALLUMPC,1);
        xlswrite('ResultsMethod3.xlsx',ALLXMPC,2);
    end
    if (SaveFileNum==4)%平方趋近于0.08
        xlswrite('ResultsMethod4.xlsx',ALLUMPC,1);
        xlswrite('ResultsMethod4.xlsx',ALLXMPC,2);
    end
end

% %% Yalmip
% %%
% Np=10;
% gamma=0.025;
% 
% 
% X_YMPC=sdpvar(3,Np+1);
% U_YMPC=sdpvar(3,Np);
% XInit=sdpvar(3,1);
% X_YMPC(:,1)=XInit;
% % assign(X_YMPC,0);
% % assign(U_YMPC,0);
% % assign(XInit,0);
% 
% Objective=0;
% Constraints=[];
% 
% Start1=tic;
% 
% for ii=1:Np
%     X_YMPC(:,ii+1)=GPRModelling(X_YMPC(:,ii),U_YMPC(:,ii),1);
%     tempVt=X_MPC(2,ii+1);
%     tempU=U_YMPC(:,ii);
%     Objective=Objective-tempVt+gamma*(tempU')*tempU;
%     Constraints=[Constraints,0<=U_YMPC(:,ii)<=2*pi];
% end
% 
% disp('定义完成');
% disp(toc(start1));
% 
% ops=sdpsettings('verbose',2);
% % optimize([Constraints,XInit == [0;0;0]],objective);
% controller=optimizer(Constraints,Objective,ops,X_YMPC(:,1),U_YMPC{:,1});
% 
% XYMPC=zeros(3,tNum);
% UYMPC=zeros(3,tNum-1);
% 
% start1=tic;
% for ii=1:tNum-1
%     
%     UYMPC(:,ii) = controller{XYMPC(:,ii)};
%     XYMPC(:,ii+1)=GPRModelling(XYMPC(:,ii),UYMPC(:,ii));
%     
% end
% disp('完成');
% disp(toc(start1));
% 
% %% fmincon
% %% 基础参数定义
% MPC.Np=10;
% MPC.delta=50/180*pi;
% MPC.N=10;
% tStep=0.01;
% tEnd=1;
% tNum=tEnd/tStep+1;
% MPC.t=0:tStep:tEnd;
% MPC.tCnt=1;
% 
% UMPC=zeros(3,tNum-1);
% XMPC=zeros(3,tNum);
% 
% XInit=[0;0;0];
% % UInit=[15/180*pi;60/180*pi;0/180*pi];
% UInit=[0;0;0];
% XMPC(:,1)=XInit;
% UMPC(:,1)=UInit;
% 
% MPC.XEveryObj=zeros(3,MPC.Np);
% MPC.XEveryObj(:,1)=XInit;
% 
% UEveryOptim=zeros(1,3*(MPC.Np-1));
% UEveryOptim(1:3)=UInit';
% 
% lb=0*ones(1,size(UEveryOptim,2));
% ub=2*pi*ones(1,size(UEveryOptim,2));
% % lb=[];
% % ub=[];
% 
% %options=optimoptions('ga','FunctionTolerance',1e-2,'ConstraintTolerance',1e-2);  
% %options=optimoptions('fmincon','Algorithm','sqp','MaxFunctionEvaluations',4000,'ConstraintTolerance',1.0000e-4,'OptimalityTolerance',1.0000e-3,'Display','notify');
% % options=optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',100,...
% %         'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6,'Display','notify');
% options=optimoptions('fmincon','Algorithm','interior-point', 'MaxFunctionEvaluations',100,...
%         'ConstraintTolerance',1e-3,'OptimalityTolerance',1e-3,'Display','iter');
%     %@NonLinearCons ,'MaxFunctionEvaluations',200,
% MPCTimeStart=tic;
% for m=1:tNum-1
%     fprintf('MPC优化：第%d个时刻\n',m);
%     [OptimU,fval,exitflag,output]=fmincon(@(UEveryOptim)ObjFunc(UEveryOptim),UEveryOptim,[],[],[],[],lb,ub,[],options);
% % [uOutPut_TotalNp,fval,exitflag]=ga(@FitFun,152,[],[],[],[],lb,ub,@NonLinearCons,options);
%     UEveryOptim=OptimU;
%     UMPC(:,m)=OptimU(1:3)';
%     XMPC(:,m+1)=GPRModelling(XMPC(:,m),UMPC(:,m),1);
%     MPC.XEveryObj(:,1)=XMPC(:,m+1);
%     UEveryOptim(1:(end-3))=UEveryOptim(4:end);
%     MPC.tCnt=MPC.tCnt+1;
% end
% 
% disp(toc(MPCTimeStart));
% 
% subplot(2,2,1)
% plot(UMPC(1,1:21));
% subplot(2,2,2)
% plot(UMPC(2,1:21));
% subplot(2,2,3)
% plot(UMPC(3,1:21));
% 
% subplot(2,2,1)
% plot(XMPC(1,1:21));
% subplot(2,2,2)
% plot(XMPC(2,1:21));
% subplot(2,2,3)
% plot(XMPC(3,1:21));