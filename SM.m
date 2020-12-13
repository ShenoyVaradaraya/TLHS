clc;
%% Symmetric motion plant
A_sm=Ap(2:5,2:5);
B_sm=Bp(2:5,2);
C_sm=Cp(2,2:5);
D_sm=Dp(2,2);
sm_sys = ss(A_sm,B_sm,C_sm,D_sm);
figure()
pzmap(sm_sys)
%% Plant Dimensions
[ns_sm,nc_sm] = size(B_sm);
no_sm = nc_sm;
%% Natural Modes: Poles (Eigenvalues), Eigenvectors
[evec,eval] = eig(A_sm)
tinit = 0;
tinc  = 0.001;
tfin  = 0.2;
t     = [tinit:tinc:tfin];   % Vector of uniformly spaced time points
u     = [0*t'];         % Set input u to zero for all time in order to generate zero input response;
% Excite Fast Instability.
% This mode 
%            is associated with a pole at s = 14.1050.
%            is associated primarily with theta_2 dot.     
x = lsim(ss(A_sm, B_sm, eye(ns_sm,ns_sm), 0*ones(ns_sm,nc_sm)), u, t, evec(:,1)); 
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
x = lsim(ss(A_sm, B_sm, eye(ns_sm,ns_sm), 0*ones(ns_sm,nc_sm)), u, t, real(evec(:,4))); 
figure; plot(t,x)
grid
title('Slow Instability')
ylabel('States')
xlabel('Time (seconds)')
pause
%% Transmission Zeros
plantzeros = tzero(ss(A_sm,B_sm,C_sm,D_sm))       % transmission zeros
%% SYSTEM TRANSFER FUNCTIONS: From u_i to x_j
Plant_zpk = zpk(ss(A_sm,B_sm,C_sm,D_sm)) % Zeros, Poles, and Gains fron u_i to x_j
%% Controllability 
rank(ctrb(A_sm,B_sm))
%% Observability 
rank(obsv(A_sm,C_sm))
%% FREQUENCY RESPONSE: Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A_sm, B_sm, C_sm, D_sm),w);
sv     = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Outputs: ;      Inputs: ')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% PLANT SVD ANALYSIS at DC
dc            =  C_sm*inv(-A_sm)*B_sm
[udc,sdc,vdc] = svd(dc)
%% Designing an LQG/LTR Compensator for symmetric motion  
[ns_sm,nc_sm] = size(B_sm);
no_sm = nc_sm;
A2 = [zeros(nc_sm) zeros(nc_sm,ns_sm); B_sm A_sm];
B2 = [eye(nc_sm);zeros(ns_sm,1)];
C2 = [zeros(nc_sm) C_sm];
D2 = [zeros(no_sm,nc_sm)];
Ll2 = inv(C_sm*inv((-A_sm))*B_sm);
Lh2 = inv(-A_sm)*B_sm*Ll2;
Llqg2 = [Ll2;Lh2];
mu = 0.1;
B_fare = (C2'*inv(mu)*C2);
sigma2 = are(A2',B_fare,(Llqg2*Llqg2'));
H2= sigma2*C2'*inv(mu);
%% Target Closed loop poles and zeros 
tclpoles2 = eig(A2-H2*C2);
tclzeros2 = tzero(A2-H2*C2,H2,C2, D2);
%% Target Open Loop Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A2,H2,C2,D2),w);
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
sv     = sigma(ss(A2-H2*C2,H2,-C2,eye(no_sm,no_sm)),w);
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
sv     = sigma(ss(A2-H2*C2,H2,C2,D2),w);
sv     = 20*log10(sv);
semilogx(w, sv)
title('Target Comp Sensitivity Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Designing Pre-Filter
z = 1.2;
fil2   = ss((tf([z],[1 z])));
s                   = j*20;
tol_s               = C2*inv(s*eye(size(A2))-A2)*H2 + D2;
[tolu, tols, tolv ] = svd(tol_s);
%% Target Closed Loop Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
[ntacl2, ntbcl2, ntccl2, ntdcl2 ] = series(fil2.a, fil2.b, fil2.c, fil2.d, A2-H2*C2,H2,C2,D2);
sv     = sigma(ss(ntacl2, ntbcl2, ntccl2, ntdcl2),w);
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
tfin  = 10.0;
t     = [tinit:tinc:tfin]';% Vector of uniformly spaced time points
r2    = 1*[ones(size(t))];
ty1   = lsim(ss(ntacl2, ntbcl2, ntccl2, ntdcl2), r2,t);

plot(t,ty1)
grid
title('Target Responses to Step Reference Command')
ylabel('(degrees)')
xlabel('Time (seconds)')
pause
%% LTR at output 
rho2 = 1e-13;
R2 = rho2*eye(1);
[G2,X2,poles2] = lqr(A2,B2,C2'*C2,R2);
%% Model Based Compensator followed by integrator bank 
Ak2 = [0*ones(nc_sm,no_sm)  G2; zeros(ns_sm+nc_sm,no_sm)  A2-B2*G2-H2*C2];
Bk2 = [zeros(nc_sm,no_sm); H2];
Ck2 = [eye(nc_sm,no_sm)  zeros(nc_sm, ns_sm+nc_sm)];
Dk2 = [zeros(no_sm,no_sm)];
kpoles = eig(Ak2)                               % Compensator Poles
kzeros = tzero(Ak2, Bk2, Ck2, Dk2)                 % Compensator Zeros

winit  = -2;
wfin   = 4;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);
sv     = sigma(ss(Ak2, Bk2, Ck2, Dk2),w);
sv     = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Form Open Loop System
[al2, bl2, cl2, dl2 ] = series(Ak2, Bk2, Ck2, Dk2, A_sm, B_sm, C_sm, D_sm);
  
olpoles = eig(al2)                          % Open Loop Poles
olzeros = tzero(al2,bl2,cl2,dl2)               % Open Loop Zeros
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(al2, bl2, cl2, dl2), w);
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
[ali2, bli2, cli2, dli2 ] = series(A_sm,B_sm,C_sm,D_sm, Ak2, Bk2, Ck2, Dk2);
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(ali2, bli2, cli2, dli2 ), w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Open Loop Singular Values at Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   
%% Form Closed Loop System
acl2    = al2-bl2*cl2;
bcl2     = bl2;
ccl2     = cl2;
dcl2     = dl2;

clpoles = eig(acl2)               % Closed Loop Poles
damp(clpoles)
clzeros = tzero(acl2,bcl2,ccl2,dcl2) % Closed Loop Zeros (r to y)
%% Sensitivity at Error
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(acl2, bcl2, -ccl2, eye(1)),w);
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
sv      = sigma(ss(ali2-bli2*cli2, bli2, -cli2, eye(1)),w);
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
sv      = sigma(ss(acl2, bcl2, ccl2, dcl2),w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values at Ouput')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Add Pre-filter to Closed Loop Transfer Function Matrix From r to y
[nacl2, nbcl2, nccl2, ndcl2 ]   = series(fil2.a, fil2.b, fil2.c, fil2.d, acl2,bcl2,ccl2,dcl2);
sv                          = sigma(ss(nacl2, nbcl2, nccl2, ndcl2),w);
sv                          = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Form Reference to Controls (r to u)
[aru2, bru2, cru2, dru2 ] = series(acl2, bcl2, -ccl2, eye(1), Ak2, Bk2, Ck2, Dk2);
[aru2, bru2, cru2, dru2 ] = series(fil2.a, fil2.b, fil2.c, fil2.d, aru2, bru2, cru2, dru2 );    
winit                 = -1;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(aru2, bru2, cru2, dru2),w);
sv                    = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Reference to Control Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Form Input Disturbance to Output
[ady2, bdy2, cdy2, ddy2 ] = series(A_sm,B_sm,C_sm,D_sm, acl2, bcl2, -ccl2, eye(1));
winit                 = -2;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(ady2, bdy2, cdy2, ddy2),w);
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
tdy2             = cdy2*inv(s*eye(size(ady2))-ady2)*bdy2 +  ddy2;
[utdy1 stdy1 vtdy1] = svd(tdy2)
%% Closed Loop Step Responses  
y1    = lsim(ss(nacl2,nbcl2,nccl2, ndcl2), r2, t); 
figure;
plot(t,y1)
grid
title('Horizontal Separation')
ylabel('Outputs \Delta x(ft)')
xlabel('Time (seconds)')
pause
u1    = lsim(ss(aru2, bru2, cru2, dru2), r2, t); 
figure; 
plot(t,u1)
grid
title('Differential Collective Control')
ylabel('Controls \Delta B_{lc} (deg)')
xlabel('Time (seconds)')
pause


