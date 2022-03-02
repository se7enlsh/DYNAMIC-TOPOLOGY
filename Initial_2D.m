%% ��ʼ2D
clc
clear

nelx = 500;                 % x����Ԫ��Ŀ
nely = 200;                  % y����Ԫ��Ŀ
sdof = 2*(nely+1)*(nelx+1); % �����ɶ�

E = 2.1e5;                  % ����ģ��
den = 7.9e-9;               % �ܶ�
nu = 0.3;                   % ���ɱ�

%% PREPARE FINITE ELEMENT ANALYSIS
% Լ��
fixeddofs = [1:(nely+1)*2];% ��߹̶�
% fixeddofs = [(nely+1)*2-1 (nely+1)*2 2*(nely+1)*(nelx+1)-1 2*(nely+1)*(nelx+1)];% ������̶�
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = E/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
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
x(1:nely,1:nelx)=1;
sK = reshape(KE(:)*x(:)',64*nelx*nely,1);
sM = reshape(ME(:)*x(:)',64*nelx*nely,1);
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
% save initial_2D Kn0 invKn0

save initial_2D Kn0 