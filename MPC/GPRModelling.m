function [XNext]=GPRModelling(X,U,ModeFlag)

%基于高斯回归的仿生蛇机器人预测函数
%
%模式0 对应GP预测（只限于数字）--->默认模式
%
%模式1 对应自己手写函数预测（支持casADi）
%
%模式2 内部函数fitrgp预测（待开发）
%
%matern 核通过函数文件内部d设置1 3 5 7
%
%输入要求：X U为3维列向量，输出X也为列向量

%选择Matern核 d=1 3 5 7

    d=1;
    
    global MPC;
    InputTest=[U',X'];

    if(nargin<3)
        ModeFlag=1;
    end

    if(ModeFlag==0)

    [ddThetaNPredAvg,ddThetaNPreCov]=gp(MPC.ddthetaNModel.hyp1,MPC.ddthetaNModel.infer,...
                                         MPC.ddthetaNModel.meanfunc,MPC.ddthetaNModel.covfunc,...
                                         MPC.ddthetaNModel.likfunc,MPC.ddthetaNModel.InputDataTrain,...
                                         MPC.ddthetaNModel.DataTrainOutput_ddthetaN,InputTest);
    [VtPreAvg,VtPreCov]=gp(MPC.VtModel.hyp2,MPC.VtModel.infer,MPC.VtModel.meanfunc,...
                                         MPC.VtModel.covfunc,MPC.VtModel.likfunc,MPC.VtModel.InputDataTrain,...
                                         MPC.VtModel.DataTrainOutput_Vt,InputTest);
    [VnPreAvg,VnPreCov]=gp(MPC.VnModel.hyp3,MPC.VnModel.infer,MPC.VnModel.meanfunc,...
                                         MPC.VnModel.covfunc,MPC.VnModel.likfunc,MPC.VnModel.InputDataTrain,...
                                         MPC.VnModel.DataTrainOutput_Vn,InputTest);
    XNext=[ddThetaNPredAvg;VtPreAvg;VnPreAvg];

    elseif ((ModeFlag==1))
        
        [ddThetaNPredAvg,ddThetaNPreCov]=GPRPredMatern(d,MPC.ddthetaNModel.hyp1.cov,...
                                 MPC.ddthetaNModel.hyp1.lik,MPC.ddthetaNModel.InputDataTrain, ...
                                 MPC.ddthetaNModel.DataTrainOutput_ddthetaN,InputTest);
        [VtPreAvg,VtPreCov]= GPRPredMatern(3,MPC.VtModel.hyp2.cov,MPC.VtModel.hyp2.lik,...
                                MPC.VtModel.InputDataTrain,MPC.VtModel.DataTrainOutput_Vt,InputTest...
                                );
        [VnPreAvg,VnPreCov]=GPRPredMatern(d,MPC.VnModel.hyp3.cov,MPC.VnModel.hyp3.lik,...
                                MPC.VnModel.InputDataTrain,MPC.VnModel.DataTrainOutput_Vn,InputTest...
                                );     
    XNext=[ddThetaNPredAvg;VtPreAvg;VnPreAvg];
    
    elseif(ModeFlag==2)
         warning('fitrgp预测模式正在开发中，请选择其他模式。');
    else
        error('GPRModelling模式选择错误，该函数未执行。');
    end

% %Fitrgp
% MPC.ddthetaNModel=load('fitrgpddthetaN.mat');
% MPC.VtModel=load('fitrgpVt.mat');
% MPC.VnModel=load('fitrgpVn.mat');
% 
% InputTest=[U',X'];
% 
% [ddThetaNPredAvg,~,~]=predict(MPC.ddthetaNModel.GPR_ddthetaN,InputTest);
% [VtPreAvg,~,~]=predict(MPC.VtModel.GPR_Vt,InputTest);
% [VnPreAvg,~,~]=predict(MPC.VnModel.GPR_Vn,InputTest);

end