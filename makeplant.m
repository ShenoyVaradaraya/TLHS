%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  makeplnt.m                                                            % 
%  Routine to generate the A,B,C,D matrices for Equal Tether Twin Lift   %
%                                                                        %
%	NOTE: Run loaddata.m prior to running this script                     %  
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Definition of variables for the Twin Lift Helicopter System           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clc;
disp('')
disp('... Calculating TLHS Data ...')
disp('')

%%
%% Mass and Moment Calculations 
M_h=W_h/g;
M_b=W_b/g;
I_b=M_b*L*L/12;
M_l=W_l/g;

%%
%% Tether Length Parameters
H_a=2*H_s*H_m/(H_s + H_m);
S=(H_m - H_s)/(H_m + H_s);
w_a_sq=g/H_a;

%%
%Normalized Quantities
H_a_hat=H_a/L;
Z_hat=Z/L;

%%
%Mass Parameters
mew=(W_l + W_b)/(2*W_h);
delta_l=W_l/(W_l +W_b);
epsl=(M_h*hp)/I_y;
ep_b=2*I_b/(M_h*L*L);

%%
%Definitons
psi=1 + ep_b + 4*mew*Z_hat*Z_hat*delta_l*(1 - delta_l);
J=inv(mew*delta_l*Z_hat*w_a_sq);
T=inv(J*psi);
D=-w_a_sq*(1 + mew + (hp + H_s)*mew*epsl + (4*T*delta_l*Z_hat/w_a_sq));
E=-(X_u + M_u*(hp + H_s));
F=(hp + H_s)*epsl*mew*w_a_sq*S*H_s - w_a_sq*H_s - 4*T*delta_l*Z_hat*H_s;
V=1 + (hp + H_s)*epsl + (1/mew) + (4*delta_l*delta_l*Z_hat*Z_hat/psi);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Definition of the Ap, Bp, Cp, and Dp Matrices used to     %
%   describe the Plant of the Twin Lift Helicopter System.    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Ap Matrix
Ap(1,1)=Z_w/(1 + mew);
Ap(2,4)=1;
Ap(3,5)=1;
Ap(4,2)=-mew*w_a_sq;
Ap(4,3)=-(g*(1 + mew) +mew*w_a_sq*hp);
Ap(4,4)=X_u;
Ap(4,6)=-2*mew*w_a_sq*S*H_s;
Ap(4,8)=-2*mew*w_a_sq*S;
Ap(5,2)=-epsl*mew*w_a_sq;
Ap(5,3)=-epsl*mew*w_a_sq*(hp + H_a);
Ap(5,4)=M_u;
Ap(5,5)=M_q;
Ap(5,6)=-2*epsl*mew*w_a_sq*S*H_s;
Ap(5,8)=-2*epsl*mew*w_a_sq*S;
Ap(6,10)=1;
Ap(7,11)=1;
Ap(8,12)=1;
Ap(9,2)=mew*w_a_sq*S/2;
Ap(9,3)=mew*w_a_sq*S*hp/2;
Ap(9,6)=-g*(1 + mew*S*H_s/H_a);
Ap(9,8)=mew*w_a_sq;
Ap(9,9)=X_u;
Ap(10,2)=mew*epsl*w_a_sq*S/2;
Ap(10,3)=mew*epsl*w_a_sq*S*hp/2;
Ap(10,6)=-epsl*mew*w_a_sq*S*H_s;
Ap(10,8)=epsl*mew*w_a_sq;
Ap(10,9)=M_u;
Ap(10,10)=M_q;
Ap(11,2)=2*T*S;
Ap(11,3)=2*T*S*hp;
Ap(11,6)=4*T*H_s;
Ap(11,7)=-4*T*H_a_hat;
Ap(11,8)=4*T;
Ap(11,11)=Z_w*T*J;
Ap(12,2)=-(mew*w_a_sq*S/2)*V;
Ap(12,3)=-(mew*w_a_sq*S/2)*V*hp;
Ap(12,6)=g*(1 + mew*S*H_s/H_a) + F;
Ap(12,7)=4*T*delta_l*Z_hat*H_a_hat;
Ap(12,8)=D;
Ap(12,9)=E;
Ap(12,10)=-M_q*(hp + H_s);
Ap(12,11)=-Z_w*T*J*delta_l*Z_hat;

%%Bp Matrix
Bp(1,1)=Z_thc/(1 + mew);
Bp(4,2)=X_blc;
Bp(5,2)=M_blc;
Bp(9,4)=X_blc;
Bp(10,4)=M_blc;
Bp(11,3)=Z_thc*T*J;
Bp(12,3)=-Z_thc*T*J*delta_l*Z_hat;
Bp(12,4)=-(X_blc + M_blc*(hp + H_s));

%%Cp Matrix
Cp(1,1)=1;
Cp(2,2)=1;
Cp(3,6)=(hp + H_s);
Cp(3,7)= Z_hat;
Cp(3,8)=1;
Cp(4,9)=1;
Cp(:,10:12)=zeros(4,3);

%%Dp Matrix
Dp=zeros(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Convert Ap, Bp, Cp, and Dp from units of radians per       %
%  second to units of degrees.                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cf=180/pi;
v1=[1,1,cf,1,cf,cf,1,1,1,cf,1,1];
con_A=diag(v1,0);
con_A_inv=inv(con_A);
v2=[cf,cf,cf,cf];
con_B=diag(v2,0);
con_B_inv=inv(con_B);
Ap=(con_A)*(Ap)*(con_A_inv);
Bp=(con_A)*(Bp)*(con_B_inv);
Cp=(Cp)*(con_A_inv);
Dp=(Dp)*(con_B_inv);

%%
%Clean Up
clear cf v1 v2 con_A con_A_inv con_B con_B_inv

%%
%% System poles and zeros
poles_sys=eig(Ap);
zeros_sys=tzero(Ap,Bp,Cp,Dp);



