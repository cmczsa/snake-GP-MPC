function [XNextAvg,XNextCov]=GPR_AvgCovModelling(X,U)

%���ڸ�˹�ع�ķ����߻�����Ԥ�⺯��
%
%ģʽ0 ��ӦGPԤ�⣨ֻ�������֣�--->Ĭ��ģʽ
%
%ģʽ1 ��Ӧ�Լ���д����Ԥ�⣨֧��casADi��
%
%ģʽ2 �ڲ�����fitrgpԤ�⣨��������
%
%matern ��ͨ�������ļ��ڲ�d����1 3 5 7
%
%����Ҫ��X UΪ3ά�����������XҲΪ������

%ѡ��Matern�� d=1 3 5 7

    d=1;
    
    global MPC;
    InputTest=[U;X]';

%     [ddThetaNPredAvg,ddThetaNPreCov]=gp(MPC.ddthetaNModel.hyp1,MPC.ddthetaNModel.infer,...
%                                          MPC.ddthetaNModel.meanfunc,MPC.ddthetaNModel.covfunc,...
%                                          MPC.ddthetaNModel.likfunc,MPC.ddthetaNModel.InputDataTrain,...
%                                          MPC.ddthetaNModel.DataTrainOutput_ddthetaN,InputTest);
%      [VtPreAvg,VtPreCov]=gp(MPC.VtModel.hyp2,MPC.VtModel.infer,MPC.VtModel.meanfunc,...
%                                          MPC.VtModel.covfunc,MPC.VtModel.likfunc,MPC.VtModel.InputDataTrain,...
%                                          MPC.VtModel.DataTrainOutput_Vt,InputTest);
%     [VnPreAvg,VnPreCov]=gp(MPC.VnModel.hyp3,MPC.VnModel.infer,MPC.VnModel.meanfunc,...
%                                          MPC.VnModel.covfunc,MPC.VnModel.likfunc,MPC.VnModel.InputDataTrain,...
%                                          MPC.VnModel.DataTrainOutput_Vn,InputTest);
%     XNextAvg=[ddThetaNPredAvg;VtPreAvg;VnPreAvg];
%     XNextCov=[ddThetaNPreCov;VtPreCov;VnPreCov];
%     
    [ddThetaNPredAvg,ddThetaNPreCov]=GPRPredMatern(d,MPC.ddthetaNModel.hyp1.cov,...
                             MPC.ddthetaNModel.hyp1.lik,MPC.ddthetaNModel.InputDataTrain, ...
                             MPC.ddthetaNModel.DataTrainOutput_ddthetaN,InputTest);
    [VtPreAvg,VtPreCov]= GPRPredMatern(d,MPC.VtModel.hyp2.cov,MPC.VtModel.hyp2.lik,...
                            MPC.VtModel.InputDataTrain,MPC.VtModel.DataTrainOutput_Vt,InputTest...
                            );
    [VnPreAvg,VnPreCov]=GPRPredMatern(d,MPC.VnModel.hyp3.cov,MPC.VnModel.hyp3.lik,...
                            MPC.VnModel.InputDataTrain,MPC.VnModel.DataTrainOutput_Vn,InputTest...
                            );     
    XNextAvg=[ddThetaNPredAvg';VtPreAvg';VnPreAvg'];
    XNextCov=[ddThetaNPreCov';VtPreCov';VnPreCov'];
    
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