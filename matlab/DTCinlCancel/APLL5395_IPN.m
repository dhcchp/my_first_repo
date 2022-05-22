% Analog CP-PLL IPN Calculation
%clc;
close all;
clear all;
%Parameter DEFINITIONS
pA = 1e-12;
nA = 1e-9;
uA = 1e-6;
mA = 1e-3;
fF = 1e-15;
pF = 1e-12;
nF = 1e-9;
uF = 1e-6;
pS = 1e-12;
MHz =1e6;
kHz = 1e3;
k = 1.38e-23;
T = 273.15 +55;

IPN_table = zeros(17,1);
PN_table = [0 0.1 0.2 0.4 0.6 0.8 1.25 1.6 1.8 3 6 10 19.3 20 30 45 80];

NTCXO = 2; %Doubler
Fxo = 50*MHz;
Fvco = 10000*MHz;
lo_div = 1;
Kvco = 2*pi*40*MHz;


%CP Noise
%f_cp = [1e3 1e4 1e5 1e6 10e6];
%Lcp_dB_scp = 1*[-196 -215.3 -220.6 -221.9 -222.3] + 10*log10(NTCXO) -10*log10(122.88/50) -3;
f_cp = [1e3 10e3 50e3 100e3 200e3 400e3 500e3 1e6 1.5e6 3e6 10e6 20e6 25e6 30e6];
Lcp_dB_scp = [-189.76 -214.18 -224.43 -226.46  -227.41 -227.80 -227.87 -227.99 -228.01 -228.02  -227.90 -227.68 -227.68 -227.68];
%XO
f_xo = [1e3 1e4 1e5 1e6 10e6];
Lref_dB = 1*[-151 -161 -170 -175 -176] + 20*log10(NTCXO);
%XO_Gen
f_xo_gen = [1e3 1e4 1e5 1e6 10e6];
Lrefgen_dB = 10*[-135 -150 -170 -173 -173] + 20*log10(NTCXO);
%VCO Noise
f_vco =   [10e3 100e3 1e6 10e6 100e6];
Lvco_dB = 1*[-79 -102 -123 -143 -163];
%Pre & NDIV Noise
f_pre = [1e4 1e5 1e6];
L_pre = 10*[-150 -170 -173] + 20*log10(NTCXO);

Icp = 5e-3;
R1 = 1.789e3;
C1 = 500*pF;
C2 = C1/10;
R3 = 2e3;
C3 =10e-12;
R4 = R3;
C4 = C3;
%R1 = 2*1/sqrt(Icp*C1/(2*pi)*Kvco/100)

%Integration BW
IPN_Fstart = 12e3;
IPN_Fstop = 20e6;

Fref = NTCXO*Fxo;
Ndiv = Fvco/Fref;
Tref = 1/Fref;

fstart = 100;
fstop = 1e8;

num_pt = 1000;
f = logspace(log10(fstart),log10(fstop),num_pt);

%__LAPACE AND Z VARIABLES
z = tf('z',Tref);
s = tf([1 0], 1, 'variable','s');
Z = exp(j*2*pi*f*Tref);
S = j*2*pi*f;

% CP Gain
Kcp = Icp/(2*pi);

Z1 = 1/(s*C3)*(R4+1/(s*C4))/(1/(s*C3)+1/(s*C4)+R4);
Zi = 1/(s*C2)*(R1+1/(s*C1))/(1/(s*C2)+1/(s*C1)+R1);
LF = Zi*(R3+Z1)/(Zi+R3+Z1) * Z1/(R3+Z1) * 1/(1+s*R4*C4);

As = Kcp*LF*Kvco/Ndiv/s;
GHs = freqs(As.num{:},As.den{:},2*pi*f);
Ts = Ndiv*GHs./(1+GHs);
Ts_mag = abs(Ts);

Hvco = 1./(1 + GHs);
Hdsm = 1e0*GHs./(1 + GHs);

Fstart = IPN_Fstart;
Fstop = IPN_Fstop;
IPN_Range = find(f >= Fstart & f <=Fstop);
%Ref Clock SSB
Lref_tot_dB = interp1(10*log10(f_xo), Lref_dB, 10*log10(f), 'linear','extrap');
So_ref_dB = Lref_tot_dB + 20*log10(abs(squeeze(Ts))) + 20*log10(1/lo_div);
IPN = 1*trapz(f(IPN_Range),10.^((So_ref_dB(IPN_Range))/10));
IPN_ref = IPN;
%Ref Buf SSB
Lrefgen_tot_dB = interp1(10*log10(f_xo_gen), Lrefgen_dB, 10*log10(f), 'linear','extrap');
So_refgen_dB = Lrefgen_tot_dB + 20*log10(abs(squeeze(Ts))) + 20*log10(1/lo_div);
IPN = 1*trapz(f(IPN_Range),10.^((So_ref_dB(IPN_Range))/10));
IPN_refgen = IPN;
%VCO SSB
L = interp1(10*log10(f_vco),Lvco_dB,10*log10(f),'linear','extrap');
So_vco_dB = L + 20*log10(abs(squeeze(Hvco))) + 20*log10(1/lo_div);
IPN = 1*trapz(f(IPN_Range),10.^((So_vco_dB(IPN_Range))/10));
IPN_vco = IPN;
%Delta Sigma Modulator SSB
Sf_dsm = 1/1*(1/(12*Fref)).*abs(((1-Z.^(-1)).^3)).^2;
L_dsm = Sf_dsm.*abs(Z.^(-1)*(2*pi)./(1-Z.^(-1))).^2;
So_dsm = L_dsm.*abs(Hdsm).^2 * (1/1)^2;
So_dsm_dB = 10*log10(So_dsm);
IPN = 1*trapz(f(IPN_Range),10.^((So_dsm_dB(IPN_Range))/10));
IPN_dsm = IPN;
%CP Noise SSB
Lcp_dB_scp_scaled = Lcp_dB_scp;% + 10*log10(Icp/(3200*uA));
L_cp_scp = interp1(10*log10(f_cp),Lcp_dB_scp_scaled, 10*log10(f),'linear','extrap');
So_cp_dB = L_cp_scp + 20*log10(abs(squeeze(Ts/Kcp))) + 20*log10(1/lo_div);
IPN_cp = 1*trapz(f(IPN_Range),10.^((So_cp_dB(IPN_Range))/10));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_C1 = (1./(S*C1));
Z_C2 = (1./(S*C2));
Z_C3 = (1./(S*C3));
Z_C4 = (1./(S*C4));

Z1_lpf = R4 + Z_C4;
Z2_lpf = (Z_C3.*(R4+Z_C4))./(Z_C3+R4+Z_C4);
Z3_lpf = R3 + Z2_lpf;
Z4_lpf =(Z_C2.*Z3_lpf)./(Z_C2+Z3_lpf);
Z5_lpf=(Z_C2.*(R1+Z_C1))./(Z_C2+R1+Z_C1);
Z6_lpf=(Z_C3.*(R3+Z5_lpf))./(Z_C3+R3+Z5_lpf);

GR1 = (Z4_lpf./(Z4_lpf+R1+Z_C1)).*(Z2_lpf./(Z2_lpf+R3)).*(Z_C4./(Z_C4 +R4));
GR2 = (Z_C4./(Z_C4+R4)).*(Z2_lpf./(Z2_lpf+R3+Z5_lpf));
GR3 = (Z_C4./(Z_C4+R4+Z6_lpf));

L_lpf = 2*k*T*(R1*(abs(GR1).^2)+R3*(abs(GR2).^2)+R4*(abs(GR3).^2));
L_lpf_dB = 10*log10(L_lpf);
So_lpf = L_lpf.*(abs((Kvco./S).*Hvco).^2);
So_lpf_dB = 10*log10(So_lpf)+20*log10(1/lo_div);
f_lpf = [1e3 10e3 50e3 100e3 500e3 1e6 1.5e6 1.6e6 1.7e6 10e6];
IPN = 1*trapz(f(IPN_Range),10.^(So_lpf_dB(IPN_Range)./10));
IPN_lpf = IPN;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IPN Plot
So_tot = 10.^(So_ref_dB./10) + 10.^(So_refgen_dB./10) + 10.^(So_vco_dB./10) + 10.^(So_dsm_dB./10) + 10.^(So_cp_dB./10) + 10.^(So_lpf_dB./10);
So_tot_dB = 10*log10(So_tot);
IPN = 1*trapz(f(IPN_Range),10.^(So_tot_dB(IPN_Range)./10));
IPN_tot = IPN;

jitter = sqrt(2*IPN)/(2*pi*Fvco)

fig_size = [200,200,910,550];
h_ipn=figure('position',fig_size);
semilogx(f,So_tot_dB,'r','LineWidth',3); grid on;
xlabel('foffset/Hz');
ylabel('PN/dBc/Hz');

axis([100 10^8 -180 -70]);

title(['PLL SSB IPN']);hold on;
semilogx(f,So_ref_dB,'k-.','LineWidth',2);
%semilogx(f,So_ref_dB+So_refgen_dB,'k-.','LineWidth',2);
semilogx(f,So_vco_dB,'m--','LineWidth',2);
semilogx(f,So_dsm_dB,'g--','LineWidth',2);
semilogx(f,So_cp_dB,'r-.','LineWidth',2);
semilogx(f,So_lpf_dB,'m-.','LineWidth',2);

legend(sprintf('tot (IPN = %0.3g dBc)',10*log10(IPN_tot)), ...
sprintf('REF (IPN = %0.3g dBc)',10*log10(IPN_ref+IPN_refgen)), ...
sprintf('VCO (IPN = %0.3g dBc)',10*log10(IPN_vco)), ...
sprintf('DSM (IPN = %0.3g dBc)',10*log10(IPN_dsm)), ...
sprintf('CP (IPN = %0.3g dBc)',10*log10(IPN_cp)), ...
sprintf('LPF (IPN = %0.3g dBc)',10*log10(IPN_lpf)), 'Location','NorthEast');

figure;margin(As);