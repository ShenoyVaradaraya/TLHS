%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  loaddata.m                                                            % 
%  Routine to load the matlab workspace with the appropriate parameters  %
%  (longitudinal dynamics)                                               % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Definition of physical parameters for the Twin Lift Helicopter System %                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Acceleration due to gravity
g=32.2;

%%
%Helicopter Parameters
W_h=14000;
I_y=5700;
hp=3.6;

%%
%Control Derivatives
X_thc=0;
X_blc=27.4;
Z_thc=340.9;
Z_blc=0;
M_thc=0;
M_blc=-47.24;

%%
%Aerodynamic derivatives at hover
X_u=-.06;
X_w=0;
X_q=0;
Z_u=0;
Z_w=-.346;
Z_q=0;
M_u=.041;
M_w=0;
M_q=-3.1;

%%
%Tether lengths
H_m=13.25;
x=1;					% tether length factor for slave
H_s=x*H_m;

%%
%%Spreader bar parameters
L=69;
W_b=644;
Z=34.5;
W_l=12000;
%
%end


