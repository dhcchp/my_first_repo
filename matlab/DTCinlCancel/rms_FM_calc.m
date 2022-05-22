clc;clear;close all;

PN =  [10e3,-99.78; 100e3,-108.73; 200e3, -112.5; 1e6,-126.02; 10e6,-147.81; 20e6,-156];

freq0 = PN(:,1);
L0 = PN(:,2);

freq = logspace(log10(freq0(1)), log10(freq0(end)), 100);
L = interp1(log10(freq0), L0, log10(freq));

semilogx(freq0,L0,'b*-');
hold on;
semilogx(freq,L,'ro-');
grid on;

L_real = (10.^(L/20));
T=L_real.^2.*freq.^2;
T2 = trapz(freq,T);
rms_FM = sqrt(2* T2 )