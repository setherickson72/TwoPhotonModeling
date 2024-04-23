close all
clear all

% Transit time calculator

m=1.44E-25; % mass in kg
T=373; % Temp in K
w=100E-6; % Beam radius
kb=1.38E-23;

v=sqrt(kb*T/m);

bw=sqrt(2*log(2))*v/w;
bw/(2*pi)