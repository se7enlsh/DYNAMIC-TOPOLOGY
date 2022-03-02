%% ��ʼ3D
clc
clear

nelx=100;
nely=30;
nelz=10;
nele = nelx*nely*nelz;               % �ܵ�Ԫ��
sdof = 3*(nelx+1)*(nely+1)*(nelz+1); % �����ɶ�

E = 2.1e5;                  % ����ģ��
den = 7.9e-9;               % �ܶ�
nu = 0.3;                   % ���ɱ�

%% PREPARE FINITE ELEMENT ANALYSIS
% Լ��
% [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);% Coordinates�����
% [iif,jf,kf] = meshgrid([0 nelx],0,0:nelz);% Coordinates�����±���
[iif,jf,kf] = meshgrid([0 nelx],0:nely,0:nelz);% Coordinates������
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
freedofs = setdiff(1:sdof,fixeddof);
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

x = ones(nely,nelx,nelz);
sK = reshape(KE(:)*x(:)',24*24*nele,1);
sM = reshape(ME(:)*x(:)',24*24*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
M = sparse(iK,jK,sM); M = (M+M')/2;
C = sparse(sdof,sdof);

%% Newmark��
dt = 0.001;                 % ʱ�䲽��
gama = 0.5;
beta = 0.25;
a0 = 1/beta/dt^2;
a1 = gama/beta/dt;
Kn = K(freedofs,freedofs)+a0*M(freedofs,freedofs)+a1*C(freedofs,freedofs); % �γ���Ч�ĸնȾ���
Kn0 = Kn;
% invKn0 = inv(Kn0);
% save initial_3D Kn0 invKn0
save initial_3D Kn0 
%% ��Ԫ�նȾ���
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
%% ��Ԫ��������
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