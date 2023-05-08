function [ddPx,ddPy,dPx,dPy,Px,Py]=SolvePcm(dtheta,ddtheta,dPx,dPy,Px,Py)
        %输入这一时刻的P，dP，输出下一时刻的P，dP，ddP
    global USR
    
%     ddX=USR.L*USR.K'*(USR.Ctheta*dtheta.^2+USR.Stheta*ddtheta)+USR.e*ddPx;
%     ddX=USR.L*USR.K'*(USR.Stheta*dtheta.^2-USR.Ctheta*ddtheta)+USR.e*ddPy;
    dX=-USR.L*USR.K'*USR.Stheta*dtheta+USR.e*dPx;%18
    dY=-USR.L*USR.K'*USR.Ctheta*dtheta+USR.e*dPy;

    Vrx=USR.Ctheta*(dX-USR.e*USR.RiverVx)+USR.Stheta*(dY-USR.e*USR.RiverVy);%39
    Vry=-USR.Stheta*(dX-USR.e*USR.RiverVx)+USR.Ctheta*(dY-USR.e*USR.RiverVy);

    fd1=-[USR.ct*USR.Ctheta,-USR.cn*USR.Stheta;USR.ct*USR.Stheta,USR.cn*USR.Ctheta]*[Vrx;Vry];%37
    fd2=-[USR.ct*USR.Ctheta,-USR.cn*USR.Stheta;USR.ct*USR.Stheta,USR.cn*USR.Ctheta]*...
            diag(sign([Vrx;Vry]))*[Vrx.^2;Vry.^2];%38

    %存下
    fd1Cell=mat2cell(fd1,[size(fd1,1)/2,size(fd1,1)/2]);
    USR.fd1x=fd1Cell{1};
    USR.fd1y=fd1Cell{2};
    fd2Cell=mat2cell(fd2,[size(fd2,1)/2,size(fd2,1)/2]);
    USR.fd2x=fd2Cell{1};
    USR.fd2y=fd2Cell{2};
    
    temp11=(USR.N*USR.m+(USR.e)'*(USR.Stheta).^2*USR.mu*USR.e);
    temp12=-(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu*USR.e;
    temp21=temp12;
    temp22=USR.N*USR.m+(USR.e)'*(USR.Ctheta).^2*USR.mu*USR.e;
    Mp=inv([temp11,temp12;temp21,temp22]);
    
    Np=[(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu,(USR.e)'*(USR.Stheta).^2*USR.mu;...
     -(USR.e)'*(USR.Ctheta).^2*USR.mu,-(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu];
    Lp=[(USR.e)'*(USR.Stheta).^2*USR.mu,-(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu;
      -(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu,(USR.e)'*(USR.Ctheta).^2*USR.mu];

    ddP=Mp*Np*[diag(dtheta),zeros(USR.N,USR.N);zeros(USR.N,USR.N),diag(dtheta)]*...
        [USR.RiverVx*USR.e;USR.RiverVy*USR.e]-Mp*Lp*[USR.L*(USR.K)'*(USR.Ctheta*...
        (dtheta).^2+USR.Stheta*ddtheta);USR.L*(USR.K)'*(USR.Stheta*(dtheta).^2....
        -USR.Ctheta*ddtheta)]+Mp*(USR.E)'*[(USR.fd1x+USR.fd2x);(USR.fd1y+USR.fd2y)];

    ddPx=ddP(1,1);
    ddPy=ddP(2,1);
    
    dPx=dPx+ddPx*USR.tStep;
    dPy=dPy+ddPy*USR.tStep;
    
    Px=Px+dPx*USR.tStep;
    Py=Py+dPy*USR.tStep;
    
%dX ddX计算方式也不一样
end