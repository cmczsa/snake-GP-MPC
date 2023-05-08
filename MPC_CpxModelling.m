function XNext=MPC_CpxModelling(X,U)
% X U 都是列向量

    global USR;
    
    persistent t;
    persistent tCnt;
    persistent theta;
    persistent dtheta;
    persistent ddtheta;
    persistent dthetaN;
    persistent Pt;
    persistent Pn;
    persistent ThistimeU;
    persistent dPy;
    persistent dPx;
    persistent Px;
    persistent Py;
    
    persistent dalpha;
    persistent dlambda;
    persistent dphi0;
    persistent PHI;
    persistent dPHI;

    if(isempty(t))
        t=0;
        tCnt=0;
        theta=zeros(USR.N,1);
        dtheta=zeros(USR.N,1);
        ddtheta=zeros(USR.N,1);
        dthetaN=0;
        Pt=0;
        Pn=0;
        Px=0;
        Py=0;
        dPx=0;
        dPy=0;
    end
    
    thetaN=X(1);
    Vt=X(2);
    Vn=X(3);
   
    alpha=U(1);
    omega=U(2);
    phi0=U(3);
    lambda=omega*t;
    
%% 
% 给前一时刻赋值
    if ~tCnt
        alphaLast=alpha;
        phi0Last=phi0;
        lambdaLast=lambda;
        dalphaLast=0;
        dlambdaLast=0;
        dphi0Last=0;
        PHILast=PHI;
        dPHILast=0;
        PxLast=0;
        PyLast=0;
    else 
        alphaLast=ThistimeU(1);
        omegaLast=ThistimeU(2);
        phi0Last=ThistimeU(3);
        lambdaLast=omegaLast.*(t-1);
        dalphaLast=dalpha;
        dlambdaLast=dlambda;
        dphi0Last=dphi0;
        PHILast=PHI;
        dPHILast=dPHI;
        PxLast=Px;
        PyLast=Py;
    end
    ThistimeU=U;



dalpha=(alpha-alphaLast)/USR.tStep;
ddalpha=(dalpha-dalphaLast)/USR.tStep;

dlambda=(lambda-lambdaLast)/USR.tStep;
% ddlambda=(dlambda-dlambdaLast)/USR.tStep;

dphi0=(phi0-phi0Last)/USR.tStep;
ddphi0=(dphi0-dphi0Last)/USR.tStep;

PHITemp=zeros(USR.N-1,1);

for jj=1:USR.N-1
    PHITemp(jj)=sin(lambda+(jj-1)*USR.delta);
end

PHI=PHITemp;
RestraintPHI=max(abs(alpha*PHI+phi0));

    if isempty(PHILast)
        PHILast=PHI;
    end 

    dPHI=(PHI-PHILast)/USR.tStep;
    ddPHI=(dPHI-dPHILast)/USR.tStep;

    USR.Stheta=diag(sin(theta));
    USR.Ctheta=diag(cos(theta));

% Pcm
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
    
    Mp_cell=mat2cell(Mp,[1 1],[1 1]);
    m11=Mp_cell{1,1};
    m12=Mp_cell{1,2};
    m21=Mp_cell{2,1};
    m22=Mp_cell{2,2};

    K1=m11*(USR.e)'*(USR.Stheta)^2*USR.mu-m12*(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu;
    K2=(-1)*m11*(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu+m12*(USR.e)'*(USR.Ctheta).^2*USR.mu;
    K3=m21*(USR.e)'*(USR.Stheta)^2*USR.mu-m22*(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu;
    K4=m21*(USR.e)'*USR.Stheta*USR.Ctheta*USR.mu-m22*(USR.e)'*(USR.Ctheta)^2*USR.mu;
    K5=(-1)*USR.m*USR.L*USR.Stheta*USR.K*USR.e-USR.L*USR.Stheta*USR.K*(USR.Stheta)^2*USR.mu*...
           USR.e-USR.L*USR.Ctheta*USR.K*USR.Stheta*USR.Ctheta*USR.mu*USR.e;
    K6=USR.m*USR.L*USR.Ctheta*USR.K*USR.e+USR.L*USR.Stheta*USR.K*USR.Stheta*USR.Ctheta*...
        USR.mu*USR.e+USR.L*USR.Ctheta*USR.K*USR.Ctheta*USR.Ctheta*USR.mu*USR.e;

    A1=(-1)*USR.L*(USR.Stheta)^2*USR.mu*(USR.K)'*USR.Stheta-USR.L*USR.Stheta*USR.Ctheta*...
        USR.mu*(USR.K)'*USR.Ctheta;
    A2=(-1)*USR.L*(USR.Stheta)^2*USR.mu*(USR.K)'*USR.Ctheta+USR.L*USR.Stheta*USR.Ctheta*...
        USR.mu*(USR.K)'*USR.Stheta;
    A4=USR.L*USR.Stheta*USR.Ctheta*USR.mu*(USR.K)'*USR.Stheta+USR.L*(USR.Ctheta)^2*USR.mu*...
        (USR.K)'*USR.Ctheta;
    A5=USR.L*USR.Stheta*USR.Ctheta*USR.mu*(USR.K)'*USR.Ctheta-USR.L*(USR.Ctheta)^2*...
        USR.mu*(USR.K)'*USR.Stheta;

    Kx=-USR.L*USR.Stheta*USR.K-K5*m11*(USR.e)'-K6*m21*(USR.e)';
    Ky=USR.L*USR.Ctheta*USR.K-K5*m12*(USR.e)'-K6*m22*(USR.e)';

    Mtheta=USR.J*USR.IN+USR.m*(USR.L)^2*USR.Stheta*USR.V*USR.Stheta+USR.m*(USR.L)^2*...
           USR.Ctheta*USR.V*USR.Ctheta-USR.L*USR.Stheta*USR.K*A1+USR.L*USR.Ctheta*USR.K*...
           A4+USR.L*K5*K1*(USR.K)'*USR.Stheta-USR.L*K5*K2*(USR.K)'*USR.Ctheta+USR.L*K6*K3*...
           (USR.K)'*USR.Stheta+USR.L*K6*K4*(USR.K)'*USR.Ctheta+USR.LAMBDA1;
    Wtheta=USR.m*(USR.L)^2*USR.Stheta*USR.V*USR.Ctheta-USR.m*(USR.L)^2*USR.Ctheta*USR.V*...
           USR.Stheta-USR.L*USR.Stheta*USR.K*A2+USR.L*USR.Ctheta*USR.K*A5+USR.L*K5*K1*...
           (USR.K)'*USR.Ctheta+USR.L*K5*K2*(USR.K)'*USR.Stheta+USR.L*K6*K3*(USR.K)'*...
           USR.Ctheta-USR.L*K6*K4*(USR.K)'*USR.Stheta;
    Vtheta=USR.LAMBDA2+USR.LAMBDA3*diag(abs(dtheta));
    Ntheta=(USR.L*USR.Stheta*USR.K*USR.Stheta*USR.Ctheta*USR.mu+USR.L*USR.Ctheta*USR.K*...
           (USR.Ctheta)^2*USR.mu-K5*K2+K6*K4)*diag(dtheta);
    Ptheta=(USR.L*USR.Stheta*USR.K*(USR.Stheta)^2*USR.mu+USR.L*USR.Ctheta*USR.K*USR.Stheta*...
            USR.Ctheta*USR.mu+K5*K1+K6*K3)*diag(dtheta);

    ddthetaN=-(USR.e'/(USR.e'*Mtheta*USR.e))*(Wtheta*dtheta.^2+Vtheta*dtheta-Ntheta*...
                      (USR.e*USR.RiverVx)-Ptheta*(USR.e*USR.RiverVy)+Kx*(USR.fd1x+USR.fd2x)...
                      +Ky*(USR.fd1y+USR.fd2y))-(USR.e'*Mtheta*USR.H*PHI/(USR.e'*Mtheta*USR.e))*...
                      ddalpha-(2*USR.e'*Mtheta*USR.H*dPHI/(USR.e'*Mtheta*USR.e))*dalpha-...
                       (USR.e'*Mtheta*USR.H*ddPHI/(USR.e'*Mtheta*USR.e))*alpha-(USR.e'*Mtheta*...
                       USR.H*USR.e_/(USR.e'*Mtheta*USR.e))*ddphi0;

    %计算下一时刻

    dthetaN=dthetaN+ddthetaN*USR.tStep;
    thetaN=thetaN+dthetaN*USR.tStep;

    ddtheta=USR.e*ddthetaN+ddalpha*USR.H*PHI+2*dalpha*USR.H*dPHI+...
                      alpha*USR.H*ddPHI+USR.H*USR.e_*ddphi0;
    dtheta=dtheta+ddtheta*USR.tStep;
    theta=theta+dtheta*USR.tStep;

    Vt=cos(thetaN)*dPx+sin(thetaN)*dPy;
    Vn=-sin(thetaN)*dPx+cos(thetaN)*dPy;
    Pt=Pt+Vt*USR.tStep;
    Pn=Pn+Vn*USR.tStep;

    XNext=[thetaN;Vt;Vn];
    
    tCnt=tCnt+1;
    t=t+USR.tStep;
end 