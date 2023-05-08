%% GPR Prepare Data 训练集和测试集时间线区分后打乱

clear
clc

% InputDataOrg=xlsread('Data.xlsx',1);
% OutputDataOrg=xlsread('Data.xlsx',2);

InputDataOrg=xlsread('Data.xlsx',1);
OutputDataOrg=xlsread('Data.xlsx',2);

% InputDataOrg=xlsread('Data1.xlsx',1);
% OutputDataOrg=xlsread('Data1.xlsx',2);

TrainDataNum=200;
EndingNum=1000;

InputDataOrg=InputDataOrg(1:EndingNum,:);
OutputDataOrg=OutputDataOrg(1:EndingNum,:);

PlotFlag1=1;

tPlot=1:size(InputDataOrg,1);
tPlot=tPlot';

thetaN=wrapToPi(InputDataOrg(:,4));
%原始数据
if(PlotFlag1)
    figure(1)
    subplot(3,2,1)
    plot(tPlot,InputDataOrg(:,1),'--','color',[0.2 0.63 0.79]);
    legend('\alpha','Fontsize',16);
    title('\bf{\it{\alpha}}','Fontsize',16);
    
    subplot(3,2,2)
    plot(tPlot,InputDataOrg(:,2),'--','color',[0.2 0.63 0.79]);
    legend('\omega','Fontsize',16)
    title('\bf{\it{\omega}}','Fontsize',16);
    subplot(3,2,3)
    plot(tPlot,InputDataOrg(:,3),'--','color',[0.2 0.63 0.79]);
    legend('\phi_0','Fontsize',16)
    title('\bf{\it{\phi_0}}','Fontsize',16);
    
    subplot(3,2,4)
    plot(tPlot,InputDataOrg(:,4),'--','color',[0.2 0.63 0.79]);
%     plot(tPlot,thetaN,'--','color',[0.2 0.63 0.79]);
    legend('$\bf{\theta_N}$','interpreter','latex')
%     title('$\bf{\ddot{\theta}_N}$','interpreter','latex','Fontsize',16);
        title('$\bf{\theta_N}$','interpreter','latex','Fontsize',16);
 
    subplot(3,2,5)
    plot(tPlot,InputDataOrg(:,5),'--','color',[0.2 0.63 0.79]);
    legend('V_t','Fontsize',16)
    title('\bf{\it{V_t}}','Fontsize',16);
 
    subplot(3,2,6)
    plot(tPlot,InputDataOrg(:,6),'--','color',[0.2 0.63 0.79]);
    legend('V_n','Fontsize',16)
    title('\bf{\it{V_n}}','Fontsize',16);
 
    figure(3)
    subplot(2,2,1)
    plot(tPlot(TrainDataNum+1:end),OutputDataOrg(TrainDataNum+1:end,1),'+');
    legend('$\bf{\ddot{\theta}_N}$','interpreter','latex','Fontsize',16)
    subplot(2,2,2)
    plot(tPlot(TrainDataNum+1:end),OutputDataOrg(TrainDataNum+1:end,2),'+');
    legend('V_t','Fontsize',16)
    subplot(2,2,3)
    plot(tPlot(TrainDataNum+1:end),OutputDataOrg(TrainDataNum+1:end,3),'+');
    legend('V_n','Fontsize',16)
    title('待预测的原始数据');
end

%% 标准化

[InputDataOrg,InputAvg,InputStd]=zscore(InputDataOrg);
[OutputDataOrg,OutputAvg,OutputStd]=zscore(OutputDataOrg);


%%存模型
% DataStandard=struct;
% DataStandard.InputAvg=InputAvg;
% DataStandard.InputStd=InputStd;
% DataStandard.OutputAvg=OutputAvg;
% DataStandard.OutputStd=OutputStd;
% save DataStandard
PlotFlag2=0;

%标准化后数据
if(PlotFlag2)
    figure(2)
    subplot(3,2,1)
    plot(tPlot,InputDataOrg(:,1),'--','color',[0.2 0.63 0.79]);
    legend('\alpha')
    subplot(3,2,2)
    plot(tPlot,InputDataOrg(:,2),'--','color',[0.2 0.63 0.79]);
    legend('\lambda')
    subplot(3,2,3)
    plot(tPlot,InputDataOrg(:,3),'--','color',[0.2 0.63 0.79]);
    legend('\phi_0')
    subplot(3,2,4)
    plot(tPlot,InputDataOrg(:,4),'--','color',[0.2 0.63 0.79]);
    legend('$\bf{\ddot{\theta}_N}$','interpreter','latex')
    subplot(3,2,5)
    plot(tPlot,InputDataOrg(:,5),'--','color',[0.2 0.63 0.79]);
    legend('V_t')
    subplot(3,2,6)
    plot(tPlot,InputDataOrg(:,6),'--','color',[0.2 0.63 0.79]);
    legend('V_n')
    
    figure(3)
    subplot(3,1,1)
    plot(tPlot,OutputDataOrg(:,1),'--','color',[0.2 0.63 0.79]);
    legend('\theta_N')
    subplot(3,1,2)
    plot(tPlot,OutputDataOrg(:,2),'--','color',[0.2 0.63 0.79]);
    legend('\V_t')
    subplot(3,1,3)
    plot(tPlot,OutputDataOrg(:,3),'--','color',[0.2 0.63 0.79]);
    legend('\V_n')
end

%% 归一化

[InputDataOrg,Input6PS]=mapminmax(InputDataOrg',-1,1);
InputDataOrg=InputDataOrg';
[OutputDataOrg,Output3PS]=mapminmax(OutputDataOrg',-4,4);
OutputDataOrg=OutputDataOrg';

figure(1)
subplot(3,2,1)
plot(tPlot,InputDataOrg(:,1),'+','color',[0.2 0.63 0.79]);
legend('\alpha')
subplot(3,2,2)
plot(tPlot,InputDataOrg(:,2),'+','color',[0.2 0.63 0.79]);
legend('\lambda')
subplot(3,2,3)
plot(tPlot,InputDataOrg(:,3),'+','color',[0.2 0.63 0.79]);
legend('\phi_0')
subplot(3,2,4)
plot(tPlot,InputDataOrg(:,4),'+','color',[0.2 0.63 0.79]);
legend('\theta_N')
subplot(3,2,5)
plot(tPlot,InputDataOrg(:,5),'+','color',[0.2 0.63 0.79]);
legend('V_t')
subplot(3,2,6)
plot(tPlot,InputDataOrg(:,6),'+','color',[0.2 0.63 0.79]);
legend('V_n')

%% 训练集和测试集分配
%分配训练集+测试集

InputDataTrainOrg=InputDataOrg(1:TrainDataNum,:);
OutputDataTrainOrg=OutputDataOrg(1:TrainDataNum,:);
InputDataTestOrg=InputDataOrg(TrainDataNum+1:end,:);
OutputDataTestOrg=OutputDataOrg(TrainDataNum+1:end,:);

InputDataTrain=zeros(size(InputDataTrainOrg,1),size(InputDataTrainOrg,2));
OutputDataTrain=zeros(size(OutputDataTrainOrg,1),size(OutputDataTrainOrg,2));
InputDataTest=zeros(size(InputDataTestOrg,1),size(InputDataTestOrg,2));
OutputDataTest=zeros(size(OutputDataTestOrg,1),size(OutputDataTestOrg,2));

% %训练集顺序
% InputDataTrain=InputDataTrainOrg;
% OutputDataTrain=OutputDataTrainOrg;
% 测试集顺序
InputDataTest=InputDataTestOrg;
OutputDataTest=OutputDataTestOrg;

%训练集乱序
RowNumTrain=randperm(size(InputDataTrainOrg,1));
for ii=1:length(InputDataTrainOrg)
    InputDataTrain(ii,:)=InputDataTrainOrg(RowNumTrain(ii),:);
    OutputDataTrain(ii,:)=OutputDataTrainOrg(RowNumTrain(ii),:);
end

%%测试集乱序
% RowNumTest=randperm(size(InputDataTestOrg,1));
% for ii=1:length(InputDataTestOrg)
%     InputDataTest(ii,:)=InputDataTestOrg(RowNumTest(ii),:);
%     OutputDataTest(ii,:)=OutputDataTestOrg(RowNumTest(ii),:);
% end

%训练集输出
DataTrainOutput_ddthetaN=OutputDataTrain(:,1);
DataTrainOutput_Vt=OutputDataTrain(:,2);
DataTrainOutput_Vn=OutputDataTrain(:,3);

%测试集输出
DataTestOutput_ddthetaN=OutputDataTest(:,1);
DataTestOutput_Vt=OutputDataTest(:,2);
DataTestOutput_Vn=OutputDataTest(:,3);

%% ----☆----☆----☆----☆----GPML Train & Predict----☆----☆----☆----☆----☆%%
%噪声标准差
NoiseSigma1=0.1;
NoiseSigma2=0.1;
NoiseSigma3=0.1;

%Start
GpTrainStart=tic;

%ddthetaN
meanfunc=[];
mean=[];
likfunc=@likGauss;
% infer=@infLaplace;
infer=@infGaussLik;%不标准化采用  标准化采用

% covfunc={'covSum', {{@covMaternard,5},{@covPPard,2},{@covRQard}}};
covfunc={@covMaternard,1};%标准化采用
% covfunc={@covMaternard,1};%不标准化采用
% cov=zeros(1,22);
cov=zeros(1,7);%标准化采用
% cov=zeros(1,2);%不标准化采用
hyp1=struct('mean',mean,'cov',cov,'lik',log(NoiseSigma1));
hyp1=minimizeGPML(hyp1,@gp,-100,infer,meanfunc,...
           covfunc,likfunc,InputDataTrain,DataTrainOutput_ddthetaN);
[ddThetaNPreAvg,ddThetaNPreCov]=gp(hyp1,infer,meanfunc,covfunc,...
                                                            likfunc,InputDataTrain,DataTrainOutput_ddthetaN,...
                                                            InputDataTest);
%存模型
save('ddthetaNModel.mat','hyp1','infer','meanfunc','covfunc','likfunc','InputDataTrain','DataTrainOutput_ddthetaN');
%Vt
meanfunc=[];%
mean=[];
likfunc=@likGauss;
% infer=@infLaplace;
infer=@infLaplace;%不标准化采用 标准化采用

% covfunc={'covSum', {{@covMaternard,3},{@covPPard,2},{@covNNone}}};
% covfunc={@covMaterniso,1};%标准化采用
% covfunc={@covMaternard,1};%不标准化采用 
covfunc={@covSEard};%标准化采用
% cov=ones(1,16);
% cov=zeros(1,2);%标准化采用
cov=zeros(1,7);%不标准化采用

hyp2=struct('mean',mean,'cov',cov,'lik',log(NoiseSigma2));
hyp2=minimizeGPML(hyp2,@gp,-100,infer,meanfunc,...
           covfunc,likfunc,InputDataTrain,DataTrainOutput_Vt);
[VtPreAvg,VtPreCov]=gp(hyp2,infer,meanfunc,covfunc,...
                                        likfunc,InputDataTrain,DataTrainOutput_Vt,InputDataTest);
%存模型
save('VtModel.mat','hyp2','infer','meanfunc','covfunc','likfunc','InputDataTrain','DataTrainOutput_Vt');
                                    
%Vn
meanfunc=[];
mean=[];
infer=@infLaplace;%标准化采用 不标准化采用
likfunc=@likGauss;
% covfunc={'covSum',{{@covMaternard,5} ,{@covPPard,2},{@covSEard}}};
% covfunc={@covMaternard,3};%标准化采用
% covfunc={@covMaternard,1};%不标准化采用
covfunc={@covSEard};%标准化采用
% cov=ones(1,21);
cov=zeros(1,7);%标准化采用
% cov=zeros(1,2);%不标准化采用

hyp3=struct('mean',mean,'cov',cov,'lik',log(NoiseSigma3));
hyp3=minimizeGPML(hyp3,@gp,-100,infer,meanfunc,...
           covfunc,likfunc,InputDataTrain,DataTrainOutput_Vn);
[VnPreAvg,VnPreCov]=gp(hyp3,infer,meanfunc,covfunc,likfunc,...
                                          InputDataTrain,DataTrainOutput_Vn,InputDataTest);                                    
%存模型
save('VnModel.mat','hyp3','infer','meanfunc','covfunc','likfunc','InputDataTrain','DataTrainOutput_Vn');


GpTrainEnd=toc(GpTrainStart);
disp('GPML Training cost:');
disp(GpTrainEnd);

%%  ----☆----☆----☆----☆----☆--GPML Reverse--☆----☆----☆----☆----☆----☆%%
%% Output 反标准化
ddThetaNPreAvg=ddThetaNPreAvg*OutputStd(1)+OutputAvg(1);
VtPreAvg=VtPreAvg*OutputStd(2)+OutputAvg(2);
VnPreAvg=VnPreAvg*OutputStd(3)+OutputAvg(3);

ddThetaNPreCov=sqrt(ddThetaNPreCov)*OutputStd(1)+OutputAvg(1);
VtPreCov=(sqrt(VtPreCov)*OutputStd(2)+OutputAvg(2));
VnPreCov=sqrt(VnPreCov)*OutputStd(3)+OutputAvg(3);

ddThetaNPreCov=ddThetaNPreCov.^2;
VtPreCov=VtPreCov.^2;
VnPreCov=VnPreCov.^2;

DataTestOutput_ddthetaN=DataTestOutput_ddthetaN*OutputStd(1)+OutputAvg(1);
DataTestOutput_Vt=DataTestOutput_Vt*OutputStd(2)+OutputAvg(2);
DataTestOutput_Vn=DataTestOutput_Vn*OutputStd(3)+OutputAvg(3);

%% Output 反归一化

OutputTestDataReverseTemp=[DataTestOutput_ddthetaN,DataTestOutput_Vt,DataTestOutput_Vn];
PredOutputReverse_AvgTemp=[ddThetaNPreAvg,VtPreAvg,VnPreAvg];
PredOutputReverse_StdTemp=[sqrt(ddThetaNPreCov),sqrt(VtPreCov),sqrt(VnPreCov)];

OutputTestDataReverse=mapminmax('reverse',OutputTestDataReverseTemp',Output3PS);
PredOutputReverse_Avg=mapminmax('reverse',PredOutputReverse_AvgTemp',Output3PS);
PredOutputReverse_Std=mapminmax('reverse',PredOutputReverse_StdTemp',Output3PS);

OutputTestDataReverse=OutputTestDataReverse';
PredOutputReverse_Avg=PredOutputReverse_Avg';
PredOutputReverse_Std=PredOutputReverse_Std';

DataTestOutput_ddthetaN=OutputTestDataReverse(:,1);
DataTestOutput_Vt=OutputTestDataReverse(:,2);
DataTestOutput_Vn=OutputTestDataReverse(:,3);

ddThetaNPreAvg=PredOutputReverse_Avg(:,1);
VtPreAvg=PredOutputReverse_Avg(:,2);
VnPreAvg=PredOutputReverse_Avg(:,3);

ddThetaNPreCov=PredOutputReverse_Std(:,1).^2;
VtPreCov=PredOutputReverse_Std(:,2).^2;
VnPreCov=PredOutputReverse_Std(:,3).^2;

%% ---☆----☆----☆----☆----☆----GPML Error----☆----☆----☆----☆----☆----☆ %%

%ddthetaN
% RealAvg=sum(abs(DataTestOutput_ddthetaN))/length(DataTestOutput_ddthetaN);
ErrorAvg=sum(abs(ddThetaNPreAvg-DataTestOutput_ddthetaN))/length(DataTestOutput_ddthetaN);
disp('ddThetaN误差(MAE)');
disp(ErrorAvg);
LowerddThetaN=ddThetaNPreAvg-2*sqrt(ddThetaNPreCov);
UpperddThetaN=ddThetaNPreAvg+2*sqrt(ddThetaNPreCov);
InRange1=DataTestOutput_ddthetaN-LowerddThetaN;
InRange2=UpperddThetaN-DataTestOutput_ddthetaN;
InRangeRate=sum((InRange1>0)&(InRange2>0))/length(DataTestOutput_ddthetaN);
disp('ThetaN在允许范围内的概率');
disp(InRangeRate);

%Vt
%RealAvg=sum(abs(DataTestOutput_Vt))/length(DataTestOutput_Vt)
ErrorAvg=sum(abs(VtPreAvg-DataTestOutput_Vt))/length(DataTestOutput_Vt);
disp('Vt误差(MAE)');
disp(ErrorAvg);
LowerVt=VtPreAvg-2*sqrt(VtPreCov);
UpperVt=VtPreAvg+2*sqrt(VtPreCov);
InRange1=DataTestOutput_Vt-LowerVt;
InRange2=UpperVt-DataTestOutput_Vt;
InRangeRate=sum((InRange1>0)&(InRange2>0))/length(DataTestOutput_Vt);
disp('Vt在允许范围内的概率');
disp(InRangeRate);

%Vn
% RealAvg=sum(abs(DataTestOutput_Vn))/length(DataTestOutput_Vn)
ErrorAvg=sum(abs(VnPreAvg-DataTestOutput_Vn))/length(DataTestOutput_Vn);
disp('Vn误差(MAE)');
disp(ErrorAvg);
LowerVn=VnPreAvg-2*sqrt(VnPreCov);
UpperVn=VnPreAvg+2*sqrt(VnPreCov);
InRange1=DataTestOutput_Vn-LowerVn;
InRange2=UpperVn-DataTestOutput_Vn;
InRangeRate=sum((InRange1>0)&(InRange2>0))/length(DataTestOutput_Vn);
disp('Vn在允许范围内的概率');
disp(InRangeRate);

%% ----☆-----☆----☆----☆----☆----GPML PLOT----☆----☆----☆----☆----☆----☆%%

SampleNum=1:length(ddThetaNPreAvg);
SampleNum=SampleNum';
%ddthetaN
figure(4)

Range95_1= [ddThetaNPreAvg+2*sqrt(ddThetaNPreCov); ...
flip(ddThetaNPreAvg-2*sqrt(ddThetaNPreCov),1)];
subplot(2,2,1)
fill([SampleNum; flip(SampleNum,1)], Range95_1, [7 7 7]/8,...%浅灰色[7 7 7]/8  浅蓝色[0.94 1 1]
    'DisplayName','95%置信区间');
hold on
plot(SampleNum,DataTestOutput_ddthetaN,'--','Color',[0.2 0.63 0.79],'Linewidth',2,...
    'DisplayName','Real Data');
plot(SampleNum,ddThetaNPreAvg,'-','Color',[1 0.5 0],'Linewidth',2,'DisplayName','Predict');
ylabel('$\bf{\it{\theta_N}}$','interpreter','latex','FontSize',16);
xlabel('测试集样本','FontSize',14);
legend('show');
title(['$\theta_N$','  ','\bf{Gaussian Process Predict}'],'interpreter','latex','Fontsize',16);

%Vt
Range95_2= [VtPreAvg+2*sqrt(VtPreCov); ...
flip(VtPreAvg-2*sqrt(VtPreCov),1)];
subplot(2,2,2)
fill([SampleNum; flip(SampleNum,1)], Range95_2, [7 7 7]/8,'DisplayName','95%置信区间');
hold on
plot(SampleNum,DataTestOutput_Vt,'--','Color',[0.2 0.63 0.79],'Linewidth',2,...
    'DisplayName','Real Data');
plot(SampleNum,VtPreAvg,'-','Color',[1 0.5 0],'Linewidth',2,'DisplayName','Predict');
ylabel('\bf{\it{V_t}}','FontSize',16);
xlabel('测试集样本','FontSize',14);
legend('show');
title(['$\bf{\it{V_t}}$','  ','\bf{Gaussian Process Predict}'],'interpreter','latex','Fontsize',16);


%Vn
Range95_3= [VnPreAvg+2*sqrt(VnPreCov); ...
flip(VnPreAvg-2*sqrt(VnPreCov),1)];
subplot(2,2,3)
fill([SampleNum; flip(SampleNum,1)], Range95_3, [7 7 7]/8,'DisplayName','95%置信区间');                                      
hold on
plot(SampleNum,DataTestOutput_Vn,'--','Color',[0.2 0.63 0.79],'Linewidth',2,...
    'DisplayName','Real Data');
plot(SampleNum,VnPreAvg,'-','Color',[1 0.5 0],'Linewidth',2,'DisplayName','Predict');
ylabel('\bf{\it{V_n}}','FontSize',16);
xlabel('测试集样本','FontSize',14);
legend('show');
title(['$\bf{\it{V_n}}$','  ','\bf{Gaussian Process Predict}'],'interpreter','latex','Fontsize',16);


 %% --------☆--------☆-------☆--------☆--------☆--------☆--------☆--------☆--%%
%% 自带Fitrgp Train +Predict+Error+Plot

GPRTraingStart=tic;
GPR_ddthetaN=fitrgp(InputDataTrain,DataTrainOutput_ddthetaN,'ResponseName','ddthetaN',...
                      'KernelFunction','matern32',...
                      'KernelParameters',[1;1],'PredictMethod','exact','FitMethod','fic',...
                     'Sigma',0.1,'Standardize',1);
GPR_Vt=fitrgp(InputDataTrain,DataTrainOutput_Vt,'ResponseName','Vt',...
                      'KernelFunction','ardmatern32', 'KernelParameters',[1;1;1;1;1;1;1],...
                     'FitMethod','exact','PredictMethod',...
                      'fic','Sigma',0.1,'Standardize',0);
GPR_Vn=fitrgp(InputDataTrain,DataTrainOutput_Vn,'ResponseName','Vn',...
                      'KernelFunction','ardmatern52', 'KernelParameters',[1;1;1;1;1;1;1],...
                      'FitMethod','exact','PredictMethod','fic','Sigma',0.1,'Standardize',0);

clc
GPRTraingEnd=toc(GPRTraingStart);
disp('Training costs;');
disp(GPRTraingEnd);

% Fitrgp Predict
[RepredTest_ddthetaN,~,LimitddThetaN]=predict(GPR_ddthetaN,InputDataTest);
LowerddThetaN=LimitddThetaN(:,1);
UpperddThetaN=LimitddThetaN(:,2);
% save('fitrgpddthetaN.mat','GPR_ddthetaN');
[RepredTest_Vt,~,LimitVt]=predict(GPR_Vt,InputDataTest);
LowerVt=LimitVt(:,1);
UpperVt=LimitVt(:,2);
% save('fitrgpVt.mat','GPR_Vt');
[RepredTest_Vn,~,LimitVn]=predict(GPR_Vn,InputDataTest);
LowerVn=LimitVn(:,1);
UpperVn=LimitVn(:,2);
% save('fitrgpVn.mat','GPR_Vn');

% fitrgp Error
ErrorPer_Pred=sum(abs(RepredTest_ddthetaN-DataTestOutput_ddthetaN))/length(DataTestOutput_ddthetaN);
ErrorPer_Real=sum(abs(DataTestOutput_ddthetaN))/length(DataTestOutput_ddthetaN);
Average=ErrorPer_Pred/ErrorPer_Real;
disp('ddThetaN误差(MAE)');
disp(ErrorPer_Pred);
InRange1=DataTestOutput_ddthetaN-LowerddThetaN;
InRange2=UpperddThetaN-DataTestOutput_ddthetaN;
InRangeRate=sum((InRange1>0)&(InRange2>0))/length(DataTestOutput_ddthetaN);
disp('ThetaN在允许范围内的概率');
disp(InRangeRate);

ErrorPer_Pred=sum(abs(RepredTest_Vt-DataTestOutput_Vt))/length(DataTestOutput_Vt);
ErrorPer_Real=sum(abs(DataTestOutput_Vt))/length(DataTestOutput_Vt);
Average=ErrorPer_Pred/ErrorPer_Real;
disp('Vt误差(MAE)');
disp(ErrorPer_Pred);
InRange1=DataTestOutput_Vt-LowerVt;
InRange2=UpperVt-DataTestOutput_Vt;
InRangeRate=sum((InRange1>0)&(InRange2>0))/length(DataTestOutput_Vt);
disp('Vt在允许范围内的概率');
disp(InRangeRate);

ErrorPer_Pred=sum(abs(RepredTest_Vn-DataTestOutput_Vn))/length(DataTestOutput_Vn);
ErrorPer_Real=sum(abs(DataTestOutput_Vn))/length(DataTestOutput_Vn);
Average=ErrorPer_Pred/ErrorPer_Real;
disp('Vn误差(MAE)');
disp(ErrorPer_Pred);
InRange1=DataTestOutput_Vn-LowerVn;
InRange2=UpperVn-DataTestOutput_Vn;
InRangeRate=sum((InRange1>0)&(InRange2>0))/length(DataTestOutput_Vn);
disp('Vn在允许范围内的概率');
disp(InRangeRate);

%Fitrgp  Plot
SampleNum=1:size(OutputDataTest,1);
figure(4)
subplot(2,2,1)
hold on 
Range1=fill([SampleNum,fliplr(SampleNum)],[LowerddThetaN',fliplr(UpperddThetaN')],...
                [0.94 1 1],'DisplayName','range');
plot(SampleNum,RepredTest_ddthetaN','b','LineWidth',2,'DisplayName','Prediction');
plot(SampleNum,DataTestOutput_ddthetaN','--','DisplayName','True Data');
% plot(SampleNum,UpperddThetaN','y--','Linewidth',1.5,'DisplayName','Upper  Limit');
% plot(SampleNum,LowerddThetaN','y--','Linewidth',1.5,'DisplayName','Lower  Limit');
ylabel('\it{\theta_N}');
xlabel('TimeStep');
legend('show');

subplot(2,2,2)
hold on 
Range2=fill([SampleNum,fliplr(SampleNum)],[UpperVt',fliplr(LowerVt')],...
                [0.94 1 1],'DisplayName','range');
plot(SampleNum,RepredTest_Vt','b','LineWidth',2,'DisplayName','Prediction');
plot(SampleNum,DataTestOutput_Vt','--','DisplayName','True Data');
% plot(SampleNum,UpperVt','y--','Linewidth',1.5,'DisplayName','Upper  Limit');
% plot(SampleNum,LowerVt','y--','Linewidth',1.5,'DisplayName','Lower  Limit');
ylabel('\it{V_t}');
xlabel('TimeStep');
legend('show');

subplot(2,2,3)
hold on 
Range2=fill([SampleNum,fliplr(SampleNum)],[UpperVn',fliplr(LowerVn')],...
                [0.94 1 1],'DisplayName','range');
plot(SampleNum,RepredTest_Vn','b','LineWidth',2,'DisplayName','Prediction');
plot(SampleNum,DataTestOutput_Vn','r+','DisplayName','True Data');
% plot(SampleNum,UpperVn','y--','Linewidth',1.5,'DisplayName','Upper  Limit');
% plot(SampleNum,LowerVn','y--','Linewidth',1.5,'DisplayName','Lower  Limit');
ylabel('\it{V_n}');
xlabel('TimeStep');
legend('show');

