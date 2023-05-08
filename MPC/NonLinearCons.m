function [c,ceq]=NonLinearCons(U)
    global MPC;
    DevMat=3*ones(1,MPC.Np-1);
    DevidedU=mat2cell(U',DevMat);
    
    TempC=zeros(MPC.Np-1,1);
    PHI=zeros(MPC.N-1,1);
   
    for ii=1:MPC.Np-1
        TempU=DevidedU{ii};
        alpha=TempU(1);
        omega=TempU(2);
        phi0=TempU(3);
        
        for jj=1:MPC.N-1
            PHI(jj)=alpha*sin(omega*MPC.t(MPC.tCnt+ii-1)+(jj-1)*MPC.delta)+phi0;
        end
        TempC(ii)=max(abs(PHI));
    end
    c=max(TempC)-90/180*pi;
    ceq=[];
end