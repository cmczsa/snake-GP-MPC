function [ymu,ys2]=GPRPredMatern(d,cov,lik,x,y,z)
%cov输入必须是行向量
%
%x为训练集输入  y为训练集输出 z为预测集
%而且三个矩阵，每行为样本，每一列为一个样本的元素
%
%仅适用于infGaussLik的收敛方式

%     if(nargin<7)
%         casFlag=0;
%     end

%     if(casFlag)
%           Ks_tempx=SX.sym('Ks_tempx',size(x,1),size(x,2));
%     end

    logsf=cov(7);
    sf=exp(logsf);

    logell=cov(1:6);%行向量

    logsn=lik;
    sn2=exp(2*logsn);

%% K
    K_mu=mean(x,1);
    
    K_muMat=repmat(K_mu,size(x,1),1);
    K_tempx=x-K_muMat;
    
%     K_tempx=x;
%     for ii=1:size(x,1)
%         K_tempx(ii,:)=x(ii,:)-K_mu;
%     end

    %计算Ax sax D2 z=[ ]
    
    ard=exp(-2*logell);
    ardMat=repmat(ard,size(K_tempx,1),1);
    K_Ax=K_tempx.*ardMat;
    
%     K_Ax=K_tempx;
%     for ii=1:size(x,1)
%         K_Ax(ii,:)=K_tempx(ii,:).*exp(-2*logell);
%     end

    K_sax = sum(K_tempx.*K_Ax,2);
    K_Az=K_Ax; 
    K_saz=K_sax;

    K_tempsaz=K_saz';
    K_temp2xAz=2*K_tempx*K_Az';% n*n

    K_tempsazMat=repmat(K_tempsaz,size(K_temp2xAz,1),1);
    temp1D2=K_tempsazMat-K_temp2xAz;
    K_saxMat=repmat(K_sax,1,size(temp1D2,2));
    temp2D2=temp1D2+K_saxMat;
    
%     temp1D2=K_temp2xAz;
%     temp2D2=temp1D2;
%     for ii=1:size(K_temp2xAz,1)
%         temp1D2(ii,:)=K_tempsaz-K_temp2xAz(ii,:);
%     end
% 
%     for ii=1:size(temp1D2,2)
%         temp2D2(:,ii)=temp1D2(:,ii)+K_sax;
%     end

    D2=max(temp2D2,0);
    %计算K0，maha文件的K

    if any(d==[1,3,5,7])
        switch d
            case 1, K0=exp(-sqrt(D2));
            case 3, K0=(1+sqrt(d*D2)).*exp(-sqrt(d*D2));
            case 5, K0=(1+sqrt(d*D2).*(1+sqrt(d*D2)/3)).*exp(-sqrt(d*D2));
            case 7, K0=(1+sqrt(d*D2).*(1+sqrt(d*D2).*(6+sqrt(d*D2))/15)).*exp(-sqrt(d*D2));
        end
    end

    %计算Scale文件的K
    sfx=sf;
    sfz=sfx;
    S= sfx*sfz';
    K = S.*K0; %n*n

    %% infGaussLik
    n=size(x,1);
    W=ones(n,1)/sn2;

    sW=sqrt(W);
    L = chol(eye(n)+sW*sW'.*K);  
    B=y.*sW;
    CholAns=L\(L'\B);%n*1
    alphaMat=CholAns.*sW;
    
    %% Kss
     D2= zeros(size(z,1),1);
    %      if d==1
    %         Kss0=exp(-sqrt(D2));%ns*1
    %      end
    if any(d==[1,3,5,7])
    switch d
        case 1, Kss0=exp(-sqrt(D2));
        case 3, Kss0=(1+sqrt(d*D2)).*exp(-sqrt(d*D2));
        case 5, Kss0=(1+sqrt(d*D2).*(1+sqrt(d*D2)/3)).*exp(-sqrt(d*D2));
        case 7, Kss0=(1+sqrt(d*D2).*(1+sqrt(d*D2).*(6+sqrt(d*D2))/15)).*exp(-sqrt(d*D2));
    end
    end
    sfx=sf;
    S=sfx.*sfx;
    sfz=sfx;
    Kss=S.*Kss0;
    
%% Ks
    ns=size(z,1);
    n=size(x,1);
    Ks_mu= (ns/(n+ns))*mean(z,1)+(n/(n+ns))*mean(x,1);
    
    Ks_ZmuMat=repmat(Ks_mu,size(z,1),1);
    Ks_tempz=z-Ks_ZmuMat;
    
%     Ks_tempz=z;
%     for ii=1:size(z,1)
%         Ks_tempz(ii,:)=z(ii,:)-Ks_mu;
%     end

    Ks_XmuMat=repmat(Ks_mu,size(x,1),1);
    Ks_tempx=x-Ks_XmuMat;
    
%     for ii=1:size(x,1)
%         Ks_tempx(ii,:)=x(ii,:)-Ks_mu;
%     end

        %计算Ax Az sax saz D2 z=样本
        
    ard=exp(-2*logell);
    KsXardMat=repmat(ard,size(Ks_tempx,1),1);
    Ks_Ax=Ks_tempx.*KsXardMat;
    
%     Ks_Ax=Ks_tempx;
%     for ii=1:size(Ks_tempx,1)
%         Ks_Ax(ii,:)=Ks_tempx(ii,:).*exp(-2*logell); %n*6
%     end

    Ks_sax = sum(Ks_tempx.*Ks_Ax,2); %n*6

    KsZardMat=repmat(ard,size(Ks_tempz,1),1);
    Ks_Az=Ks_tempz.*KsZardMat;
    
%     Ks_Az=Ks_tempz;
%     for ii=1:size(z,1)
%         Ks_Az(ii,:)=Ks_tempz(ii,:).*exp(-2*logell);%ns*6
%     end

    Ks_saz = sum(Ks_tempz.*Ks_Az,2); %ns*6

    Ks_tempsaz=Ks_saz';
    Ks_temp2xAz=2*Ks_tempx*Ks_Az';

    Ks_tempsazMat=repmat(Ks_tempsaz,size(Ks_temp2xAz,1),1);
    temp3D2=Ks_tempsazMat-Ks_temp2xAz;
    Ks_saxMat=repmat(Ks_sax,1,size(temp3D2,2));
    temp4D2=temp3D2+Ks_saxMat;

%     temp3D2=Ks_temp2xAz;
%     temp4D2=temp3D2;
%     for ii=1:size(Ks_temp2xAz,1)
%         temp3D2(ii,:)=Ks_tempsaz-Ks_temp2xAz(ii,:);
%     end
% 
%     for ii=1:size(temp3D2,2)
%         temp4D2(:,ii)=temp3D2(:,ii)+Ks_sax;
%     end

    D2=max(temp4D2,0);%n*ns

        %计算K0，maha文件的K
    %     if d==1
    %         Ks0=exp(-sqrt(D2));
    %     end
     if any(d==[1,3,5,7])
        switch d
            case 1, Ks0=exp(-sqrt(D2));
            case 3, Ks0=(1+sqrt(d*D2)).*exp(-sqrt(d*D2));
            case 5, Ks0=(1+sqrt(d*D2).*(1+sqrt(d*D2)/3)).*exp(-sqrt(d*D2));
            case 7, Ks0=(1+sqrt(d*D2).*(1+sqrt(d*D2).*(6+sqrt(d*D2))/15)).*exp(-sqrt(d*D2));
        end
     end    

    %计算Scale文件的K
    sfx=sf;
    sfz=sf;
    S= sfx*sfz';
    Ks= S.*Ks0; 
    
%% 计算输出
    N = size(alphaMat,2);  
    Fmu=Ks'*full(alphaMat);
    %fmu=sum(Fmu,2)/N;

    V  = L'\(repmat(sW,1,ns).*Ks);
    
%     for ii=1:size(V,1)
%             V(ii)=V(ii)^2;
%     end
%     fs2=Kss-sum(V.*V,1)';
    
    fs2= Kss - sum(V.*V,1)';   
    fs2= max(fs2,0); 
    Fs2 = repmat(fs2,1,N);

    Ymu=Fmu;
    Ys2=Fs2+sn2;
    ymu=sum(Ymu,2)/N; 
    ys2 = sum(Ys2,2)/N;   
end