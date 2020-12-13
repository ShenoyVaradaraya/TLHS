clc;
A_asm=Ap(6:12,6:12);
B_asm=Bp(6:12,3:4);
C_asm=Cp(3:4,6:12);
D_asm=Dp(3:4,3:4);
poles_asm=eig(A_asm);
zeros_asm=tzero(A_asm,B_asm,C_asm,D_asm);
asm_sys = ss(A_asm,B_asm,C_asm,D_asm);
figure()
pzmap(asm_sys)
%% Plant Dimensions
[ns_asm,nc_asm] = size(B_asm);
no_asm = nc_asm;
%% Natural Modes: Poles (Eigenvalues), Eigenvectors
[evec,eval] = eig(A_asm)
tinit = 0;
tinc  = 0.001;
tfin  = 0.2;
t     = [tinit:tinc:tfin];   % Vector of uniformly spaced time points
u     = [0*t' 0*t'];         % Set input u to zero for all time in order to generate zero input response;
% Excite Fast Instability.
% This mode 
%            is associated with a pole at s = 14.1050.
%            is associated primarily with theta_2 dot.     
x = lsim(ss(A_asm, B_asm, eye(ns_asm,ns_asm), 0*ones(ns_asm,nc_asm)), u, t, evec(:,1)); 
figure; plot(t,x)
grid
title('Fast Instability')
ylabel('States')
xlabel('Time (seconds)')
pause
%
% Excite Slow Instability.
% This mode 
%            is associated with a poles at s = 4.5299.
%            is associated primarily with theta_1 and theta_1 dot
%
x = lsim(ss(A_asm, B_asm, eye(ns_asm,ns_asm), 0*ones(ns_asm,nc_asm)), u, t, real(evec(:,4))); 
figure; plot(t,x)
grid
title('Slow Instability')
ylabel('States')
xlabel('Time (seconds)')
pause
%% Transmission Zeros
plantzeros = tzero(ss(A_asm,B_asm,C_asm,D_asm))       % transmission zeros
%% SYSTEM TRANSFER FUNCTIONS: From u_i to x_j
Plant_zpk = zpk(ss(A_asm,B_asm,C_asm,D_asm)) % Zeros, Poles, and Gains fron u_i to x_j
%% Controllability 
rank(ctrb(A_asm,B_asm))
%% Observability 
rank(obsv(A_asm,C_asm))
%% FREQUENCY RESPONSE: Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A_asm, B_asm, C_asm, D_asm),w);
sv     = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Outputs: ;      Inputs: ')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% PLANT SVD ANALYSIS at DC
dc            =  C_asm*inv(-A_asm)*B_asm
[udc,sdc,vdc] = svd(dc)
%% Designing an LQG/LTR Compensator for anti-symmetric motion
A3 = [zeros(nc_asm,ns_asm+nc_asm); B_asm A_asm];
B3 = [eye(nc_asm);zeros(ns_asm,nc_asm)];
C3 = [zeros(nc_asm) C_asm];
D3 = [zeros(no_asm,nc_asm)];
Ll3 = inv(C_asm*inv((-A_asm))*B_asm);
Lh3 = inv(-A_asm)*B_asm*Ll3;
Llqg3 = [Ll3;Lh3];
mu = 1.5;
B_fare3 = (C3'*inv(mu*eye(2))*C3);
sigma3 = are(A3',B_fare3,(Llqg3*Llqg3'));
H3= sigma3*C3'*inv(mu*eye(2));
%% Target Closed loop poles and zeros 
tclpoles3 = eig(A3-H3*C3);
tclzeros3 = tzero(A3-H3*C3,H3,C3, D3);
%% Target Open Loop Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A3,H3,C3,D3),w);
tlsv   = 20*log10(sv);
semilogx(w, tlsv)
title('Target Open Loop Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Target Sensitivity Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A3-H3*C3,H3,-C3,eye(no_asm,no_asm)),w);
sv     = 20*log10(sv);
semilogx(w, sv)
title('Target Sensitivity Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Target Complementary Sensitivity Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A3-H3*C3,H3,C3,D3),w);
sv     = 20*log10(sv);
semilogx(w, sv)
title('Target Comp Sensitivity Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Designing Pre-Filter
z = 1.2;
fil3   =  ss(tf({z 0; 0 z}, {[1 z] 1; 1 [1 z]}));
s                   = j*20;
tol_s               = C3*inv(s*eye(size(A3))-A3)*H3 + D3;
[tolu, tols, tolv ] = svd(tol_s);
%% Target Closed Loop Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
[ntacl3, ntbcl3, ntccl3, ntdcl3 ] = series(fil3.a, fil3.b, fil3.c, fil3.d, A3-H3*C3,H3,C3,D3);
sv     = sigma(ss(ntacl3, ntbcl3, ntccl3, ntdcl3),w);
sv     = 20*log10(sv);
semilogx(w, sv)
title('Target Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('T_{ry}, Singular Values (dB)')
pause
%% Target Step Responses
tinit = 0;
tinc  = 0.005;
tfin  = 100.0;
t     = [tinit:tinc:tfin]';% Vector of uniformly spaced time points
r1    = -1*[ones(size(t)) zeros(size(t))];
r2    = 5*[zeros(size(t)) ones(size(t))];
ty1   = lsim(ss(ntacl3, ntbcl3, ntccl3, ntdcl3), r1,t);
ty2   = lsim(ss(ntacl3, ntbcl3, ntccl3, ntdcl3), r2,t);

plot(t,ty1(:,1),'b', t,ty1(:,2),'r')
grid
title('Target Responses to Step Reference Command')
ylabel('')
xlabel('Time (seconds)')
pause


plot(t,ty2(:,1),'b', t,ty2(:,2),'r')
grid
title('Target Responses to Step Reference Command')
ylabel('')
xlabel('Time (seconds)')
pause
%% LTR at output 
rho3 = 0.001;
R3 = rho3*eye(nc_asm);
[G3,X3,poles3] = lqr(A3,B3,C3'*C3,R3);
%% Model Based Compensator followed by integrator bank 
Ak3 = [0*ones(nc_asm,no_asm)  G3; zeros(ns_asm+nc_asm,no_asm)  A3-B3*G3-H3*C3];
Bk3 = [zeros(nc_asm,no_asm); H3];
Ck3 = [eye(nc_asm,no_asm)  zeros(nc_asm, ns_asm+nc_asm)];
Dk3 = [zeros(no_asm,no_asm)];
kpoles = eig(Ak3)                               % Compensator Poles
kzeros = tzero(Ak3, Bk3, Ck3, Dk3)                 % Compensator Zeros

winit  = -2;
wfin   = 4;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);
sv     = sigma(ss(Ak3, Bk3, Ck3, Dk3),w);
sv     = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Form Open Loop System
[al3, bl3, cl3, dl3 ] = series(Ak3, Bk3, Ck3, Dk3, A_asm, B_asm, C_asm, D_asm);
  
olpoles = eig(al3)                          % Open Loop Poles
olzeros = tzero(al3,bl3,cl3,dl3)               % Open Loop Zeros
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(al3, bl3, cl3, dl3), w);
sv      = 20*log10(sv);
figure; semilogx(w, sv, w, tlsv)
%clear sv
title('Open Loop Singular Values at Error (Recovered and Target)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   


figure; semilogx(w, sv)
%clear sv
title('Open Loop Singular Values at Error')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
[ali3, bli3, cli3, dli3 ] = series(A_asm,B_asm,C_asm,D_asm, Ak3, Bk3, Ck3, Dk3);
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(ali3, bli3, cli3, dli3 ), w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Open Loop Singular Values at Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   
%% Form Closed Loop System
acl3    = al3-bl3*cl3;
bcl3     = bl3;
ccl3     = cl3;
dcl3     = dl3;

clpoles = eig(acl3)               % Closed Loop Poles
damp(clpoles)
clzeros = tzero(acl3,bcl3,ccl3,dcl3) % Closed Loop Zeros (r to y)
%% Sensitivity at Error
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(acl3, bcl3, -ccl3, eye(2)),w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values at Error')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Sensitivity at Input
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(ali3-bli3*cli3, bli3, -cli3, eye(2)),w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values at Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Complementary Sensitivity
winit   = -1;
wfin    = 2;
nwpts   = 200;
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(acl3, bcl3, ccl3, dcl3),w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values at Ouput')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Add Pre-filter to Closed Loop Transfer Function Matrix From r to y
[nacl3, nbcl3, nccl3, ndcl3 ]   = series(fil3.a, fil3.b, fil3.c, fil3.d, acl3,bcl3,ccl3,dcl3);
sv                          = sigma(ss(nacl3, nbcl3, nccl3, ndcl3),w);
sv                          = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Form Reference to Controls (r to u)
[aru3, bru3, cru3, dru3 ] = series(acl3, bcl3, -ccl3, eye(2), Ak3, Bk3, Ck3, Dk3);
[aru3, bru3, cru3, dru3 ] = series(fil3.a, fil3.b, fil3.c, fil3.d, aru3, bru3, cru3, dru3 );    
winit                 = -1;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(aru3, bru3, cru3, dru3),w);
sv                    = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Reference to Control Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Form Input Disturbance to Output
[ady3, bdy3, cdy3, ddy3 ] = series(A_asm,B_asm,C_asm,D_asm, acl3, bcl3, -ccl3, eye(2));
winit                 = -2;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(ady3, bdy3, cdy3, ddy3),w);
sv                    = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Input Disturbance to Output Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Input Disturbance Analysis at DC
s = j*7
tdy3            = cdy3*inv(s*eye(size(ady3))-ady3)*bdy3 +  ddy2;
[utdy1 stdy1 vtdy1] = svd(tdy3)
%% Closed Loop Step Responses  
y1    = lsim(ss(nacl3,nbcl3, nccl3, ndcl3), r1, t); 
figure;
plot(t,y1(:,1))
grid
title('Response to Step Reference Command for Load Deviation from Center ')
ylabel('Outputs Load Deviation from Center , x_{L}-\Sigma x')
xlabel('Time (seconds)')
pause
%return

y2    = lsim(ss(nacl3,nbcl3, nccl3, ndcl3), r2, t); 
figure;
plot(t,y2(:,2))
grid
title('Response to Step Reference Command for Average horizontal velocity')
ylabel('Outputs \Sigma \dot x')
xlabel('Time (seconds)')
pause
%return


u1    = lsim(ss(aru3, bru3, cru3, dru3), r1, t); 
figure; plot(t,u1(:,2))
grid
title('Response to Step Reference Command for Differential Cyclic Control')
ylabel('Controls \theta_{c}')
xlabel('Time (seconds)')
pause
%return


u2    = lsim(ss(aru3, bru3, cru3, dru3), r2, t); 
figure;plot(t,u2(:,2))
grid
title('Response to Step Reference Command for Average Cyclic Control')
ylabel('Controls \Sigma \beta_{lc}')
xlabel('Time (seconds)')
pause
