function [ymu,ys2]=GPRPredRQ(cov,lik,x,y,z)
%cov输入必须是行向量
%仅适用于infGaussLik
%预测变量z 必须是行向量ns*6

%     if(nargin<6)
%         casFlag=0;
%     end
    
%     if(casFlag)
%          Ks_tempx=SX.sym('Ks_tempx',size(x,1),size(x,2));
%     end
    
    logsf=cov(7);
%     logsf=4.1304;
    sf=exp(logsf);
    
    logell=cov(1:6);
%     logell=[-0.0158, 0.8299,0.2337,4.5145,0.7350,-1.0882];

    logalpha=cov(8);
%     logalpha=1.0013;
    alpha=exp(logalpha);
    
    logsn=lik;
    sn2=exp(2*logsn);
%     sn2=exp(2* 3.9715);
    
 %% z=[ ]下D2 K
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
%     for ii=1:size(temp1D2,2)
%         temp2D2(:,ii)=temp1D2(:,ii)+K_sax;
%     end

    D2=max(temp2D2,0);
    %计算K0，maha文件的K
    K0=(1+0.5*D2/alpha).^(-alpha);
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
    Kss0=(1+0.5*D2/alpha).^(-alpha);%ns*1
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
    

%         Ks_Ax=Ks_tempx;
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
%     for ii=1:size(temp3D2,2)
%         temp4D2(:,ii)=temp3D2(:,ii)+Ks_sax;
%     end

    D2=max(temp4D2,0);%n*ns
    
    %计算K0，maha文件的K
    Ks0=(1+0.5*D2/alpha).^(-alpha);
    %计算Scale文件的K
    sfx=sf;
    sfz=sf;
    S= sfx*sfz';
    Ks= S.*Ks0; 
    
%% 计算输出
    % if issparse(alphaMat)
    %     nz=alphaMat~=0;
    % else
    %     nz=true(size(alphaMat,1),1);
    % end
    N = size(alphaMat,2);  
    Fmu=Ks'*full(alphaMat);
    fmu=sum(Fmu,2)/N;

    V  = L'\(repmat(sW,1,ns).*Ks);
%     fs2= Kss - sum(V.*V,1)';   
    fs2= Kss - sum(V.*V,1)'; 
    fs2= max(fs2,0); 
    Fs2 = repmat(fs2,1,N);

    Ymu=Fmu;
    Ys2=Fs2+sn2;
    ymu=sum(Ymu,2)/N; 
    ys2 = sum(Ys2,2)/N; 
end 