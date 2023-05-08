function [Vax,Vay]=Rotation(theta,VWater_Global)
% VWater_Global，theta都是列向量 
    VWater_Linki=zeros(length(VWater_Global),length(theta));
    for ii=1:length(theta)
        VWater_Linki(:,ii)=[cos(theta(ii)),-sin(theta(ii));sin(theta(ii)),cos(theta(ii))]*VWater_Global;
    end
%     Va=[diag(VWater_Linki(1,:));diag(VWater_Linki(2,:))];
    Vax=diag(VWater_Linki(1,:));
    Vay=diag(VWater_Linki(2,:));
end