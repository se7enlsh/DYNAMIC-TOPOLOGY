%% 动态拓扑2D(ESL POD)
clc
clear

nelx = 500;                 % x方向单元数目
nely = 200;                  % y方向单元数目
sdof = 2*(nely+1)*(nelx+1); % 总自由度
volfrac = 0.5;              % 体积分数
penal = 3;                  % 惩罚权重
rmin = 2;                   % 过滤半径
ft = 1;                     % 1-敏度过滤 2-密度过滤
Maxloop = 400;              % 最大迭代步数

a = 0.5;                    % 单元密度值变化量
b = 0.01;                   % 变化单元所占百分比

E0 = 2.1e5;                 % 弹性模量
Emin = E0*1e-9;
den = 7.9e-9;               % 密度
nu = 0.3;                   % 泊松比

dt = 0.001;                 % 时间步长
t = 0.2;                    % 总时间
nt = floor(t/dt);           % 时间步数

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
ME = (den/9)*[4 0 2 0 1 0 2 0
              0 4 0 2 0 1 0 2
              2 0 4 0 2 0 1 0
              0 2 0 4 0 2 0 1
              1 0 2 0 4 0 2 0
              0 1 0 2 0 4 0 2
              2 0 1 0 2 0 4 0
              0 2 0 1 0 2 0 4];
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% 约束
% fixeddofs = [1:(nely+1)*2];% 左边固定
fixeddofs = [(nely+1)*2-1 (nely+1)*2 2*(nely+1)*(nelx+1)-1 2*(nely+1)*(nelx+1)];% 三点弯固定
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% freedofs = [1:2*(nely+1)*(nelx+1)];
% 载荷
F = sparse(sdof,nt+1);
t1 = 0:dt:t;
% F(2*(nely+1)*(nelx/2)+2,1:1:nt+1) = 100*cos(50*t1);% 上边中间
% F(2*(nely+1)*(nelx/2+1),1:1:nt+1) = -100*cos(50*t1);% 下边中间
% F(2*(nely+1)*(nelx+1),1:1:nt+1) = -100*sin(50*t1);% 右下端点
% loaddofs = [1:2:(nely+1)*2];% 左边
% F(loaddofs,1:1:nt+1) = repmat(-100,[length(loaddofs),nt+1]);% 左边
% loaddofs = [2*(nely+1)*(nelx+1)-1:-2:2*(nely+1)*nelx+1];% 右边
% F(loaddofs,1:1:nt+1) = repmat(-100*sin(50*t1),[length(loaddofs),1]);% 右边

% % 随机激励
% % Frand = -150+300*rand(1,nt+1);
% Frand = -300+600*rand(1,nt+1);
% % % save FRAND Frand
% % % load FRAND
% smooth = smooth(Frand,10,'rlowess');% 'moving''lowess''loess''sgolay''rlowess''rloess'
% save suiji1 smooth
load suiji
% % plot(Frand);
% % hold on
plot(t1,smooth);
F(2*(nely+1)*(nelx/2+1),1:1:nt+1) = smooth;
% F(2*(nely+1)*(nelx+1),1:1:nt+1) = smooth;
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

x0(1:nely,1:nelx)=0;
x1(1:nely,1:nelx)=1;
x(1:nely,1:nelx)=volfrac;

loop1 = 0;
while length(find(abs(x1-x0)>=a)) > b*nelx*nely
    tic
    loop1 = loop1+1;
    x0=x1;
    sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    sM = reshape(ME(:)*(Emin+x(:)'.^penal),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    M = sparse(iK,jK,sM); M = (M+M')/2;
%     % 加质量点
%     massdofs = [2*(nely+1)*(nelx+1)-1:-1:2*(nely+1)*nelx+1];
%     M(massdofs,massdofs) = M(massdofs,massdofs) + 10;
    
    C = sparse(sdof,sdof);

    d = sparse(sdof,nt+1);   % 位移
    vel = sparse(sdof,nt+1); % 速度
    
%     %加初速度
%     veldofs = [2*(nely+1)*(nelx+1)-1:-2:2*(nely+1)*nelx+1];
%     vel(veldofs,1) = repmat(10000000,[nely+1,1]);
    
    acc = sparse(sdof,nt+1); % 加速度
    
    [d,vel,acc] = Newmark(K,M,C,F,d,vel,acc,dt,nt,freedofs);
    toc
    time = 0:dt:t;
    plot(time,d(2*(nely+1)*(nelx/2+1),:))
    xlabel('time(seconds)')
    ylabel('disp')
    %saveas(gcf,[num2str(loop1),'th iter ','DOF',num2str(2*(nely+1)*(nelx/2+1)),'.bmp'],'bmp');
    %% 等效静载
    
%     % 完整等效载荷
%     Feq = K*d;

%     % 方法1 选择峰值点
%     Z = d(2*(nely+1)*(nelx/2+1),:);
%     it = 20;
%     k = 1;
%     for i=1:it:nt+1
%         Q(1,k)=Z(i);
%         Q(2,k)=i;
%         k=k+1;
%     end
%     k=1;
%     for i=1:size(Q,2)-2
%         while (Q(1,i+1)>=Q(1,i)&&Q(1,i+1)>=Q(1,i+2))||(Q(1,i+1)<=Q(1,i)&&Q(1,i+1)<=Q(1,i+2))
%             EQ(1,k)=Q(2,i+1);
%             k=k+1;
%             break
%         end
%     end
%     Feq = K*d(:,EQ);
%     
    % 方法2 POD
     
    Feq = K*d;
    [U,S,~] = svds(Feq);
    plot(diag(S));
    xlabel ( 'Proper Orthogonal Mode' ) ;
    ylabel ( 'Singular value' ) ;
    set(gca,'xtick',1:1:6);
    saveas(gcf,'奇异值随机2维下中.bmp');
    nPOD = 2; % 选择的本证正交基数目 
    Feq = U(:,1:nPOD);
  
    n = size(Feq,2); % 等效静载数目

    %% START ITERATION
    xPhys = x;
    loop = 0;
    change = 1;
    U = zeros(sdof,n); 
    
    while change > 0.01 && loop < Maxloop
        tic
        loop = loop + 1;
        %% FE-ANALYSIS
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs,:) = K(freedofs,freedofs)\Feq(freedofs,:);
        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        c = 0;
        dc = 0;
        for i = 1:n
            Ui = U(:,i);
            ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
            c = c + sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
            dc = dc - penal*(E0-Emin)*xPhys.^(penal-1).*ce;
        end
            dv = ones(nely,nelx);
        %% FILTERING/MODIFICATION OF SENSITIVITIES
        if ft == 1
            dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
        elseif ft == 2
            dc(:) = H*(dc(:)./Hs);
            dv(:) = H*(dv(:)./Hs);
        end
        %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e9; move = 0.2;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
            if ft == 1
                xPhys = xnew;
            elseif ft == 2
                xPhys(:) = (H*xnew(:))./Hs;
            end
            if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
        end
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        %% PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
            mean(xPhys(:)),change);
        %% PLOT DENSITIES
        toc
        colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    end    
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    saveas(gcf,[num2str(loop1),'.bmp'],'bmp');
    % 过滤
    e1=0.5;
    x1=x;
    for ely=1:nely
        for elx=1:nelx
            if x1(ely,elx)<=e1
                x1(ely,elx)=0;
            else
                x1(ely,elx)=1;
            end
        end
    end
    
end
%% Newmark法，f为freedofs
function [disp,vel,acc]=Newmark(K,M,C,F,disp,vel,acc,dt,nt,f)  
gama = 0.5;
beta = 0.25;
a0 = 1/beta/dt^2;
a1 = gama/beta/dt;
a2 = 1/beta/dt;
a3 = 1/2/beta-1;
a4 = gama/beta-1;
a5 = dt/2*(gama/beta-2);
a6 = dt*(1-gama);
a7 = gama*dt;
K1 = K(f,f)+a0*M(f,f)+a1*C(f,f); % 形成有效的刚度矩阵
%invK1 = inv(K1);
%[L,U] = lu(K1);
acc(f,1) = M(f,f)\(F(f,1)-K(f,f)*disp(f,1)-C(f,f)*vel(f,1)); % t=0时刻的加速度
for i = 2:nt+1
    f1 = F(f,i)+M(f,f)*(a0*disp(f,i-1)+a2*vel(f,i-1)+a3*acc(f,i-1))...
        +C(f,f)*(a1*disp(f,i-1)+a4*vel(f,i-1)+a5*acc(f,i-1)); % 有效载荷
    %disp(f,i) = invK1*f1;
    
%     UD = L\f1;
%     disp(f,i) = U\UD;
    
    disp(f,i) = K1\f1;    
    
    acc(f,i) = a0*(disp(f,i)-disp(f,i-1))-a2*vel(f,i-1)-a3*acc(f,i-1);
    vel(f,i) = vel(f,i-1)+a6*acc(f,i-1)+a7*acc(f,i);
end
end

