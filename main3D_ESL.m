%% 动态拓扑3D(ESL POD)
clc
clear

nelx = 4;
nely = 1;
nelz = 2;
volfrac = 0.5;              % 体积分数
penal = 3;                  % 惩罚权重
rmin = 2;                   % 过滤半径
% 迭代控制参数
a = 0.5;                    % 单元密度值变化量
b = 0.01;                   % 变化单元所占百分比
maxloop = 200;              % 最大迭代步数
tolx = 0.015;               % 拓扑迭代终止准则
displayflag = 1;            % 可视化：1-每步显示；0-只显示最后一步
% 材料属性
E0 = 2.1e5;                 % 弹性模量
Emin = E0*1e-9;             % 空材料弹性模量
den = 7.9e-9;               % 密度
nu = 0.3;                   % 泊松比
% Newmark积分参数
dt = 0.001;                 % 时间步长
t = 0.2;                    % 总时间
nt = floor(t/dt);           % 时间步数

%% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;               % 总单元数
sdof = 3*(nelx+1)*(nely+1)*(nelz+1); % 总自由度
KE = lk_H8(nu);
ME = lm_H8(den);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% 约束
% [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);% Coordinates左端面
% [iif,jf,kf] = meshgrid([0 nelx],0,0:nelz);% Coordinates两端下边线
[iif,jf,kf] = meshgrid([0 nelx],0:nely,0:nelz);% Coordinates两端面
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
freedofs = setdiff(1:sdof,fixeddof);
% 载荷
% [il,jl,kl] = meshgrid(nelx, 0, 0:nelz);% Coordinates右端下边线
% [il,jl,kl] = meshgrid(nelx/2, nely, 0:nelz);% Coordinates上中线
[il,jl,kl] = meshgrid(nelx/2, 0, 0:nelz);% Coordinates下中线
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % DOFs
F = sparse(sdof,nt+1);
t1 = 0:dt:t;
% F(loaddof,1:1:nt+1) = repmat(-10000*sin(50*t1),[length(loaddof),1]);
load suiji
F(loaddof,1:1:nt+1) = repmat(smooth',[length(loaddof),1]);

%% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);

x0 = zeros(nely,nelx,nelz);
x1 = ones(nely,nelx,nelz);
x = repmat(volfrac,[nely,nelx,nelz]);
loop1 = 0;

while length(find(abs(x1-x0)>=a)) > b*nelx*nely
    loop1 = loop1+1;
    x0=x1;
    sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),24*24*nele,1);
    sM = reshape(ME(:)*(Emin+x(:)'.^penal),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    M = sparse(iK,jK,sM); M = (M+M')/2;
    C = sparse(sdof,sdof);

    d = sparse(sdof,nt+1);   % 位移
    vel = sparse(sdof,nt+1); % 速度
    acc = sparse(sdof,nt+1); % 加速度
    tic
    d = Newmark(K,M,C,F,d,vel,acc,dt,nt,freedofs);
    toc
%     d0 = d;
%     save d100x25x20 d0
    time = 0:dt:t;
%     plot(time,d(3*(nely+1)*(nelx+1)*(nelz+1)-1,:))
    plot(time,d(3*(nely+1)*(nelx/2)*(nelz/2)-1,:))
    xlabel('time(seconds)')
    ylabel('disp')
    %saveas(gcf,[num2str(loop1),'th iter ','DOF',num2str(2*(nely+1)*(nelx/2+1)),'.bmp'],'bmp');
    %% 等效静载
    
%     % 完整等效载荷
%     Feq = K*d;

%     % 方法1 选择峰值点
% %     Z = d(sdof,:);
%     seldof = 3*(nely+1)*(nelx/2)*(nelz/2)-1;
%     Z = d(seldof,:);    
%     it = 5;
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
    
    % 方法2 POD
    Feq = K*d;
    [U,S,~] = svds(Feq);
    plot(diag(S));
    xlabel ( 'Proper Orthogonal Mode' ) ;
    ylabel ( 'Singular value' ) ;
    set(gca,'xtick',1:1:6);
    saveas(gcf,'奇异值随机3维.bmp');
    nPOD = 2; % 选择的本证正交基数目 
    Feq = U(:,1:nPOD);
    
    n = size(Feq,2); % 等效静载数目
    %% INITIALIZE ITERATION
    
    xPhys = x;
    loop = 0;
    change = 1;
    U = zeros(sdof,n);
    % START ITERATION
    
    while change > tolx && loop < maxloop
        tic
        loop = loop+1;
        % FE-ANALYSIS
        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        U(freedofs,:) = K(freedofs,freedofs)\Feq(freedofs,:);
        % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        c = 0;
        dc = 0;
        for i = 1:n
            Ui = U(:,i);
            ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),[nely,nelx,nelz]);
            c = c + sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
            dc = dc - penal*(E0-Emin)*xPhys.^(penal-1).*ce;
        end
        dv = ones(nely,nelx,nelz);
        % FILTERING AND MODIFICATION OF SENSITIVITIES
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
        % OPTIMALITY CRITERIA UPDATE
        l1 = 0; l2 = 1e9; move = 0.2;
        while (l2-l1)/(l1+l2) > 1e-3
            lmid = 0.5*(l2+l1);
            xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
            xPhys(:) = (H*xnew(:))./Hs;
            if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
        end
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        % PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
        % PLOT DENSITIES
        toc
        if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
    end
    clf; display_3D(xPhys);
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


%% 单元刚度矩阵
function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
%% 单元质量矩阵
function [ME] = lm_H8(den)
ME = (den/27)*[ 8 0 0 4 0 0 2 0 0 4 0 0 4 0 0 2 0 0 1 0 0 2 0 0;
                0 8 0 0 4 0 0 2 0 0 4 0 0 4 0 0 2 0 0 1 0 0 2 0;
                0 0 8 0 0 4 0 0 2 0 0 4 0 0 4 0 0 2 0 0 1 0 0 2;
                4 0 0 8 0 0 4 0 0 2 0 0 2 0 0 4 0 0 2 0 0 1 0 0;
                0 4 0 0 8 0 0 4 0 0 2 0 0 2 0 0 4 0 0 2 0 0 1 0;
                0 0 4 0 0 8 0 0 4 0 0 2 0 0 2 0 0 4 0 0 2 0 0 1;
                2 0 0 4 0 0 8 0 0 4 0 0 1 0 0 2 0 0 4 0 0 2 0 0;
                0 2 0 0 4 0 0 8 0 0 4 0 0 1 0 0 2 0 0 4 0 0 2 0;
                0 0 2 0 0 4 0 0 8 0 0 4 0 0 1 0 0 2 0 0 4 0 0 2;
                4 0 0 2 0 0 4 0 0 8 0 0 2 0 0 1 0 0 2 0 0 4 0 0;
                0 4 0 0 2 0 0 4 0 0 8 0 0 2 0 0 1 0 0 2 0 0 4 0;
                0 0 4 0 0 2 0 0 4 0 0 8 0 0 2 0 0 1 0 0 2 0 0 4;
                4 0 0 2 0 0 1 0 0 2 0 0 8 0 0 4 0 0 2 0 0 4 0 0;
                0 4 0 0 2 0 0 1 0 0 2 0 0 8 0 0 4 0 0 2 0 0 4 0;
                0 0 4 0 0 2 0 0 1 0 0 2 0 0 8 0 0 4 0 0 2 0 0 4;
                2 0 0 4 0 0 2 0 0 1 0 0 4 0 0 8 0 0 4 0 0 2 0 0;
                0 2 0 0 4 0 0 2 0 0 1 0 0 4 0 0 8 0 0 4 0 0 2 0;
                0 0 2 0 0 4 0 0 2 0 0 1 0 0 4 0 0 8 0 0 4 0 0 2;
                1 0 0 2 0 0 4 0 0 2 0 0 2 0 0 4 0 0 8 0 0 4 0 0;
                0 1 0 0 2 0 0 4 0 0 2 0 0 2 0 0 4 0 0 8 0 0 4 0;
                0 0 1 0 0 2 0 0 4 0 0 2 0 0 2 0 0 4 0 0 8 0 0 4;
                2 0 0 1 0 0 2 0 0 4 0 0 4 0 0 2 0 0 4 0 0 8 0 0;
                0 2 0 0 1 0 0 2 0 0 4 0 0 4 0 0 2 0 0 4 0 0 8 0;
                0 0 2 0 0 1 0 0 2 0 0 4 0 0 4 0 0 2 0 0 4 0 0 8 ];
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
    t1 = cputime;
%     UD = L\f1;
%     disp(f,i) = U\UD;
    disp(f,i) = K1\f1;
    t2 = cputime;
    t = t2 -t1;
    acc(f,i) = a0*(disp(f,i)-disp(f,i-1))-a2*vel(f,i-1)-a3*acc(f,i-1);
    vel(f,i) = vel(f,i-1)+a6*acc(f,i-1)+a7*acc(f,i);
end
end
%% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;            % User-defined unit element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
set(gcf,'Name','ISO display','NumberTitle','off');
for k = 1:nelz
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            if (rho(j,i,k) > 0.5)  % User-defined display density threshold
                vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end
