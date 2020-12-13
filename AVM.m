clc;
%% Average vertical motion plant
A_avm=Ap(1,1);
B_avm=Bp(1,1);
C_avm=Cp(1,1);
D_avm=Dp(1,1);
poles_avm=eig(A_avm);
zeros_avm=tzero(A_avm,B_avm,C_avm,D_avm);
avm_sys = ss(A_avm,B_avm,C_avm,D_avm);
figure()
pzmap(avm_sys)
%% Plant Dimensions
[ns_avm,nc_avm] = size(B_avm);
no_avm = nc_avm;
%% Eigenvectors and Eigenvalues
[evec,eval] = eig(A_avm);
%% Modal Analysis 
%
tinit = 0;
tinc  = 0.001;
tfin  = 0.2;
t     = [tinit:tinc:tfin];   % Vector of uniformly spaced time points
u     = [0*t'];         % Set input u to zero for all time in order to generate zero input response;
                             % i.e. response to an initial condition x_o.
                             
x = lsim(ss(A_avm, B_avm, eye(ns_avm,ns_avm), 0*ones(ns_avm,nc_avm)), u, t, evec); 
figure; plot(t,x)
grid
title('Fast Instability')
ylabel('States')
xlabel('Time (seconds)')
pause
% Excite Slow Instability.
x = lsim(ss(A_avm,B_avm , eye(ns_avm,ns_avm), 0*ones(ns_avm,nc_avm)), u, t, real(evec)); 
figure; plot(t,x)
grid
title('Slow Instability:')
ylabel('States (deg, deg/sec)')
xlabel('Time (seconds)')
%% Transmission Zeros
plantzeros = tzero(ss(A_avm,B_avm,C_avm,D_avm))       % transmission zeros
%% Plant Transfer Funtion from u to x
Plant_zpk = zpk(ss(A_avm,B_avm,C_avm,D_avm));
%% Controllability 
rank(ctrb(A_avm,B_avm))
%% Observability 
rank(obsv(A_avm,C_avm))
%% Frequency Response 
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A_avm, B_avm, C_avm, D_avm),w);
sv     = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Outputs: \Sigma dz(ft/sec) ; Inputs: \Sigma \theta_{c}')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Pland SVD analysis at DC
dc            =  C_avm*inv(-A_avm)*B_avm
[udc,sdc,vdc] = svd(dc)

%% Designing an LQG/LTR Compensator for average vertical motion 
A1 = [zeros(nc_avm,nc_avm+ns_avm);B_avm A_avm];
B1 = [eye(nc_avm); zeros(nc_avm)];
C1 = [zeros(nc_avm) C_avm];
D1 = [zeros(no_avm,nc_avm)];
Ll1 = inv(C_avm*inv((-A_avm))*B_avm);
Lh1 = inv(-A_avm)*B_avm*Ll1;
Llqg1 = [Ll1;Lh1];
mu = 0.1;
sigma1 = are(A1',(C1'*inv(mu)*C1),(Llqg1*Llqg1'));
H1 = sigma1*C1'*inv(mu);
%% Target Closed loop poles and zeros 
tclpoles = eig(A1-H1*C1);
tclzeros = tzero(A1-H1*C1,H1,C1, D1);
%% Target Open Loop Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
sv     = sigma(ss(A1,H1,C1,D1),w);
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
sv     = sigma(ss(A1-H1*C1,H1,-C1,eye(no_avm,no_avm)),w);
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
sv     = sigma(ss(A1-H1*C1,H1,C1,D1),w);
sv     = 20*log10(sv);
semilogx(w, sv)
title('Target Comp Sensitivity Singular Values at Output')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Designing Pre-Filter
z = 1.2;
fil1    = ss((tf([z],[1 z])));
s                   = j*20;
tol_s               = C1*inv(s*eye(size(A1))-A1)*H1 + D1;
[tolu, tols, tolv ] = svd(tol_s);
%% Target Closed Loop Singular Values
winit  = -1;
wfin   =  2;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);   % Form vector of logarithmically spaced freq points
[ntacl1, ntbcl1, ntccl1, ntdcl1 ] = series(fil1.a, fil1.b, fil1.c, fil1.d, A1-H1*C1,H1,C1,D1);
sv     = sigma(ss(ntacl1, ntbcl1, ntccl1, ntdcl1),w);
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
tfin  = 4.0;
t     = [tinit:tinc:tfin]';                   % Vector of uniformly spaced time points
r1    = 5*[ones(size(t))];
ty1   = lsim(ss(ntacl1, ntbcl1, ntccl1, ntdcl1), r1,t);

plot(t,ty1)
grid
title('Target Responses to Step Reference Command')
ylabel('(degrees)')
xlabel('Time (seconds)')
pause
%% LTR at output 
rho1 = 1e-13;
R1 = rho1*eye(1);
[G1,X1,poles1] = lqr(A1,B1,C1'*C1,R1);
%% Model Based Compensator followed by integrator bank 
Ak1 = [0*ones(nc_avm,no_avm)  G1; zeros(ns_avm+nc_avm,no_avm)  A1-B1*G1-H1*C1];
Bk1 = [zeros(nc_avm,no_avm); H1];
Ck1 = [eye(nc_avm,no_avm)  zeros(nc_avm, ns_avm+nc_avm)];
Dk1 = [zeros(no_avm,no_avm)];
kpoles = eig(Ak1)                               % Compensator Poles
kzeros = tzero(Ak1, Bk1, Ck1, Dk1)                 % Compensator Zeros

winit  = -2;
wfin   = 4;
nwpts  = 200;  
w      = logspace(winit,wfin,nwpts);
sv     = sigma(ss(Ak1, Bk1, Ck1, Dk1),w);
sv     = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause
%% Form Open Loop System
[al1, bl1, cl1, dl1 ] = series(Ak1, Bk1, Ck1, Dk1, A_avm, B_avm, C_avm, D_avm);
  
olpoles = eig(al1)                          % Open Loop Poles
olzeros = tzero(al1,bl1,cl1,dl1)               % Open Loop Zeros
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(al1, bl1, cl1, dl1), w);
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
[ali1, bli1, cli1, dli1 ] = series(A_avm, B_avm, C_avm, D_avm, Ak1, Bk1, Ck1, Dk1);
    
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(ali1, bli1, cli1, dli1 ), w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Open Loop Singular Values at Input')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   
%% Form Closed Loop System
acl1     = al1-bl1*cl1;
bcl1     = bl1;
ccl1     = cl1;
dcl1     = dl1;

clpoles = eig(acl1)               % Closed Loop Poles
damp(clpoles)
clzeros = tzero(acl1,bcl1,ccl1,dcl1) % Closed Loop Zeros (r to y)
%% Sensitivity at Error
winit   = -1;
wfin    = 2;
nwpts   = 200;  
w       = logspace(winit,wfin,nwpts);
sv      = sigma(ss(acl1, bcl1, -ccl1, eye(1)),w);
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
sv      = sigma(ss(ali1-bli1*cli1, bli1, -cli1, eye(1)),w);
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
sv      = sigma(ss(acl1, bcl1, ccl1, dcl1),w);
sv      = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values at Ouput')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Add Pre-filter to Closed Loop Transfer Function Matrix From r to y
[nacl1, nbcl1, nccl1, ndcl1 ]   = series(fil1.a, fil1.b, fil1.c, fil1.d, acl1,bcl1,ccl1,dcl1);
sv                          = sigma(ss(nacl1, nbcl1, nccl1, ndcl1),w);
sv                          = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Closed Loop Singular Values (r to y)')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Form Reference to Controls (r to u)
[aru1, bru1, cru1, dru1 ] = series(acl1, bcl1, -ccl1, eye(1), Ak1, Bk1, Ck1, Dk1);
[aru1, bru1, cru1, dru1 ] = series(fil1.a, fil1.b, fil1.c, fil1.d, aru1, bru1, cru1, dru1 );    
winit                 = -1;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(aru1, bru1, cru1, dru1),w);
sv                    = 20*log10(sv);
figure; semilogx(w, sv)
%clear sv
title('Reference to Control Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 
%% Form Input Disturbance to Output
[ady1, bdy1, cdy1, ddy1 ] = series(A_avm, B_avm, C_avm, D_avm, acl1, bcl1, -ccl1, eye(1));
winit                 = -2;
wfin                  = 2;
nwpts                 = 200;
w                     = logspace(winit,wfin,nwpts);
sv                    = sigma(ss(ady1, bdy1, cdy1, ddy1),w);
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
tdy1              = cdy1*inv(s*eye(size(ady1))-ady1)*bdy1 +  ddy1;
[utdy1 stdy1 vtdy1] = svd(tdy1)
%% Closed Loop Step Responses     
y1    = lsim(ss(nacl1,nbcl1, nccl1, ndcl1), r1, t); 
figure; plot(t,y1, t, ty1)
plot(t,y1)
grid
title('Average Vertical Velocity')
ylabel('Outputs \Sigma dz(ft/sec)')
xlabel('Time (seconds)')
pause
u1    = lsim(ss(aru1, bru1, cru1, dru1), r1, t); 
figure; 
plot(t,u1)
grid
title('Average Collective Control')
ylabel('Controls \Sigma \theta_{c}')
xlabel('Time (seconds)')
pause


