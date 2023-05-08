function Obj=ObjFunc(U)

    global MPC;
    DevMat=3*ones(1,MPC.Np-1);
    DevidedU=mat2cell(U',DevMat);
%     disp('ObjFun´ÎÊý£º')
%     disp(MPC.tCnt);

    for ii=1:MPC.Np-1
        MPC.XEveryObj(:,ii+1)=GPRModelling(MPC.XEveryObj(:,ii),DevidedU{ii},1);
    end
    
    Obj=-sum(MPC.XEveryObj(2,:));
    
end