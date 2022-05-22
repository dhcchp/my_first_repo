%% General Description
% This is the TCAS-I DTC based ADPLL architecture.
% The DTC is variable slope DTC with fix pattern INL. The Main pattern is
% to be cancelled.



% Version 3, Peng CHEN, 2021.12.12

%% Initial Parameter Settings
clc;clear;close all;
Sim_Time = 2^19; % The number is in FREF cycles.
settle_sample = 20e6; % The unit is in CKV cycles.

%% Defining parameters
fref = 50e6;
Tref = 1/fref;
f0 = 2.4e9;
T0 = 1/f0;
fcw = 50+1/2^6+1*1/3;%When FCW_F is too small, the INL calibration does not work. It is related to the PLL LF coefficients.
FCW_I = floor(fcw); FCW_F = mod(fcw,1);
fv = fcw*fref;
Tv = 1/fv;
delta_TDC = 50e-12;
N_TDC_avg = 128;
KDCO = 80e3;
Mdcorun = 2;
alpha = 2^(-6);
rho = 2^(-14);
%%% DCO 1/f^2 noise
f_offset = 1e6;
S_DCO_offset_dB = -110;
S_DCO_offset = 10^(S_DCO_offset_dB/10);
sigma_w = f_offset/fv*sqrt(Tv)*sqrt(S_DCO_offset);
delta_Tw = sigma_w*randn(round(Sim_Time*fcw*Mdcorun),1);
%%% DCO thermal noise
S_DCO_thermal_dB = -160;
S_DCO_thermal = 10^(S_DCO_thermal_dB/10);
sigma_j = 1/(2*pi*fv)*sqrt(S_DCO_thermal*fv);
delta_Tj = sigma_j*randn(round(Sim_Time*fcw*Mdcorun),1);
%%% Reference noise
S_ref_dB = -170;
S_ref = 10^(S_ref_dB/10);
sigma_ref = 1/(2*pi*fref)*sqrt(S_ref*fref);
delta_Tref = sigma_ref*randn(Sim_Time,1);

%% Variables initialization
% Variables changing at fref rate
Rr = zeros(round(Sim_Time),1);
otw = zeros(size(Rr));
tR = zeros(size(Rr));
tR_norm = zeros(size(Rr));
Delta_TDEV = zeros(size(Rr));
trise = zeros(size(Rr));
Tv_current = zeros(size(Rr));
Tv_est = zeros(size(Rr));
epsilon = zeros(size(Rr));
phie = zeros(size(Rr));
mash_out = 0*order2_mash(0.5,Sim_Time);
DCW = zeros(size(Rr));
DTC_delay = zeros(size(Rr));
% Variables changing at fv rate
Rv = zeros(round(Sim_Time*fcw*Mdcorun),1);
TDEV = zeros(size(Rv));
TDEVpn = zeros(size(Rv));
tckv = zeros(size(Rv));
tckv_norm = zeros(size(Rv));
tckv_div = zeros(size(Rv));
tckv_period = zeros(size(Rv));
% number of bits for the key parameters;Here it is the frac bits 
NB_dco =  1*5; % We need <= 10KHz DCO resolution; But this value reduces when the BW reduces.
step_DCO = 1;
I_path = 0;
yiir1 = 0; yiir2 = 0; yiir3 = 0; yiir4 = 0;lambda1 = 1/2^1;lambda2 = 0.8;lambda3=0.75;lambda4=1/2^0;
tmp = zeros(size(Rr));
tmp2 = zeros(size(Rr));
tmp3 = zeros(size(Rr));
tmp4 = zeros(size(Rr));
tmp5 = zeros(size(Rr));
PHRF = 0;
DTC_reso = 0.5e-12;
DTC_gain = Tv/DTC_reso;
DTC_gain_est = 0.97*DTC_gain;
rng(17);
DTC_DNL = 1e-13*randn(1,floor(1.7*floor(Tv/DTC_reso)));
%DTC_INL = cumsum(DTC_DNL);
DTC_INL = 10e-12*sin(1.5*pi*[1:length(DTC_DNL)]/length(DTC_DNL));
%DTC_INL = 10e-12*( -([1:length(DTC_DNL)]-0.7*length(DTC_DNL)/2).^2+length(DTC_DNL)^2/4)/length(DTC_DNL)^2*4;
EN_correction = 0;
prbs = 10*randi(2,1,length(Rr));
tckv_adjust = 0; tckv_adjust_prev = 0;
DTC_gain2_est = 0;
DTC_gain3_est = 0;
DTC_gain4_est = 0;
DTC_gain5_est = 0;

%% Main Loop
if step_DCO == 1
    TDEV(1) = 0 + Delta_TDEV(1);
    TDEVpn(1) = 0 + delta_Tj(1) - 0 + delta_Tw(1);
    tckv(1) = 1*T0 - TDEV(1) + TDEVpn(1);
    tckv_norm(1) = 1*Tv;
    tckv_div(1) = tckv(1) - step_DCO*Tv;
    tckv_period(1) =  tckv(1);
end
for step = 1:Sim_Time
    PHRF = PHRF + FCW_F;
    OV_PHRF = floor(PHRF);
    PHRF = mod(PHRF,1);
    if ((prbs(step)==2)&& (PHRF<0.5))
        DCW(step) = floor((1-PHRF-0.5)*DTC_gain_est + 1*(1-PHRF-0.5)^2*DTC_gain2_est + 1*(1-PHRF-1*0.5)^3*DTC_gain3_est )+1 ;
        DCW_residue = (1-PHRF-0.5)*DTC_gain_est+1  - DCW(step) ;
    elseif ((prbs(step)==2) && (PHRF>= 0.5))
        DCW(step) = floor((1-PHRF+0.5)*DTC_gain_est + 1*(1-PHRF-0.5)^2*DTC_gain2_est + 1*(1-PHRF-1*0.5)^3*DTC_gain3_est )+1 ;
        DCW_residue = (1-PHRF+0.5)*DTC_gain_est+1  - DCW(step);
    else
        DCW(step) = floor((1-PHRF)*DTC_gain_est + 1*(1-PHRF-0.5)^2*DTC_gain2_est + 1*(1-PHRF-1*0.5)^3*DTC_gain3_est + 1*(1-PHRF-0.5)^4*DTC_gain4_est + 1*(1-PHRF-0.5)^5*DTC_gain5_est)+1 ;
        DCW_residue = (1-PHRF)*DTC_gain_est+1  - DCW(step);
    end
    
    DTC_delay(step) = DCW(step)*DTC_reso+DTC_INL(DCW(step));
    
    if step == 1  
        Delta_TDEV(1) = T0 - 1/(f0 + otw(1)*KDCO);        
        tR(1) = step * Tref + delta_Tref(1)+DTC_delay(1);        
        tR_norm(step) = step * Tref;        
    else        
        Delta_TDEV(step) = T0 - 1/(f0 + otw(step)*KDCO);        
        tR(step) = step * Tref + delta_Tref(step)+DTC_delay(step);
        tR_norm(step) = step * Tref;
    end
    CNT_cnt = 0;
    
    tckv_adjust_prev = tckv_adjust;
    if ((prbs(step)==2)&& (PHRF<0.5))
        tckv_adjust = -0.5;
    elseif ((prbs(step)==2) && (PHRF>= 0.5))
        tckv_adjust = +0.5;
    else
        tckv_adjust = 0;
    end    
    
    while ((tckv(step_DCO)+tckv_adjust*tckv_period(step_DCO)) < tR(step))
        step_DCO = step_DCO + 1;
        TDEV(step_DCO) = TDEV(step_DCO - 1) + Delta_TDEV(step);
        TDEVpn(step_DCO) = TDEVpn(step_DCO-1) + delta_Tj(step_DCO)-delta_Tj(step_DCO-1) + delta_Tw(step_DCO);
        tckv(step_DCO) = step_DCO*T0 -(TDEV(step_DCO)) + TDEVpn(step_DCO);
        tckv_norm(step_DCO) = step_DCO*Tv;
        tckv_div(step_DCO) = tckv(step_DCO) - step_DCO*Tv;
        tckv_period(step_DCO) =  tckv(step_DCO) - tckv(step_DCO - 1);
        CNT_cnt = CNT_cnt + 1;
    end
    trise(step) = (-tR(step) + tckv(step_DCO)+1*tckv_adjust*tckv_period(step_DCO) - 100e-12);
    Tv_est = 1.0*Tv; 
    
    epsilon(step) = (floor(trise(step)/delta_TDC)-0.5)*delta_TDC/Tv_est - (CNT_cnt - FCW_I - OV_PHRF) - 0*DCW_residue/DTC_gain_est ...
                    ;
    phie(step) =  epsilon(step);
    
    DTC_gain_est = DTC_gain_est + 2*1e-0*phie(step)*(1-PHRF-0.5);
    if (DTC_gain_est <= 0.5*DTC_gain) DTC_gain_est = 0.5*DTC_gain; elseif (DTC_gain_est >= 1.5*DTC_gain) DTC_gain_est = 1.5*DTC_gain; end 
    %DTC_gain_est = 783.52;
    
    DTC_gain2_est = DTC_gain2_est + 6*1e-0*phie(step)*((1-PHRF-0.5)^2-0.^2);
    if (DTC_gain2_est <= -100) DTC_gain2_est = -100; elseif DTC_gain2_est >= 100 DTC_gain2_est = 100; end
    %DTC_gain2_est = 64.47;

    DTC_gain3_est = DTC_gain3_est + 26*1e-0*phie(step)*(1-PHRF-0.5)^3;
    if (DTC_gain3_est <= -50) DTC_gain3_est = -50; elseif DTC_gain3_est >= 50 DTC_gain3_est = 50; end
    %DTC_gain3_est = 22.43;
    
    DTC_gain4_est = DTC_gain4_est + 100*1e-0*phie(step)*(1-PHRF-0.5)^4;
    if (DTC_gain4_est <= -100) DTC_gain4_est = -100; elseif DTC_gain4_est >= 100 DTC_gain4_est = 100; end
    
    DTC_gain5_est = DTC_gain5_est + 100*1e-0*phie(step)*(1-PHRF-0.5)^5;
    if (DTC_gain5_est <= -100) DTC_gain5_est = -100; elseif DTC_gain5_est >= 100 DTC_gain5_est = 100; end
    
    tmp(step) = DTC_gain_est;
    tmp2(step) = DTC_gain2_est;
    tmp3(step) = DTC_gain3_est;
    tmp4(step) = DTC_gain4_est;
    tmp5(step) = DTC_gain5_est;
    
    yiir1 = (1-lambda1)*yiir1 + lambda1*(phie(step));
    yiir2 = (1-lambda2)*yiir2 + lambda2*yiir1;
    yiir3 = (1-lambda3)*yiir3 + lambda3*yiir2;
    yiir4 = (1-lambda4)*yiir4 + lambda4*yiir3;
    P_path = alpha*yiir4;
    I_path = I_path + phie(step);
    LFout = P_path+rho*I_path;
    
    otw(step+1) = LFout * fref/KDCO;
    otw(step+1) = floor(otw(step+1)*2^NB_dco)/2^NB_dco;
end

%% Plot the Inst frequency deviation of CKV
if 1 % Determine plot or not.
    fckv_div = 1./tckv_period(1:settle_sample) - fv;
    fckv_div_avg = zeros(size(fckv_div));
    N_fckv_div_avg = 25;
    for num = 1: settle_sample
        if num < N_fckv_div_avg
            fckv_div_avg(num) = mean(fckv_div(1:num+N_fckv_div_avg));
        elseif (num > (settle_sample-N_fckv_div_avg))
            fckv_div_avg(num) = mean(fckv_div(num-N_fckv_div_avg+1:settle_sample));
        else
            fckv_div_avg(num) = mean(fckv_div(num-N_fckv_div_avg+1:num+N_fckv_div_avg));
        end
    end
    Fig = figure('Color',[1 1 1]);
    plot(1e6*tckv(1:settle_sample),fckv_div/1e6);
    hold on
    plot(1e6*tckv(1:settle_sample),fckv_div_avg/1e6,'r-');
    grid on;
    xlabel('Time (us)','fontsize',16);
    ylabel('Freq Deviation (MHz)','fontsize',16);
    xlim([0*Tv*1e6 settle_sample*Tv*1e6])
    set(gca,'YTick',-100:5:100);%This is just to tune the y-axis step, rather than the range.
    set(gca,'fontsize',16)
    set(findobj(gca,'Type','line'),'LineWidth',2)
    set(gcf,'color','w');
end

%% Plot the Inst phase error phie
if 1
    Fig = figure('Color',[1 1 1]);
    %phie = (tckv_div(1:Sim_Time)-mean(tckv_div(1:Sim_Time)))./Tv*2*pi;
    plot(tR_norm(1:Sim_Time)*1e6,phie);
    hold on;grid on;
    xlabel('Time (us)','fontsize',16)
    ylabel('Phase Deviation \phi_e','fontsize',16);
    xlim([0*Tv settle_sample*Tv*1e6])
    set(gca,'YTick',0:0.5:5);
    set(gca,'fontsize',16)
    set(findobj(gca,'Type','line'),'LineWidth',2)
    set(gcf,'color','w');
end

%% Plot the Inst Period deviation of CKV
if 1
    Fig = figure('Color',[1 1 1]);
    period_div = (tckv_div(1:step_DCO)-mean(tckv_div(settle_sample:step_DCO)));
    plot(tckv_norm(1:step_DCO)*1e6,period_div*1e12);
    hold on;grid on;
    xlabel('Time (us)','fontsize',20)
    ylabel('Period Deviation (ps)','fontsize',20);
    xlim([0*Tv settle_sample*Tv*1e6]);
    set(gca,'YTick',-2500:500:2500);
    set(gca,'fontsize',16);
    set(findobj(gca,'Type','line'),'LineWidth',2);
    set(gcf,'color','w');
end

%% S-domain model
figure;
rbw = 1e3;fstep = rbw;
f = 0:rbw:fv-rbw;
PN_tdc = (2*pi)^2/12*(delta_TDC/Tv)^2/fref * ones(size(f));
PN_vco = S_DCO_offset * (f_offset ./f).^2;
%Hol = (alpha + rho*fref./(j*2*pi*f)).*fref./(j*2*pi*f);
Hol = (alpha + rho * fref./(1 - exp(-2*pi*j*f/fref)) /fref ) .* fref ./ (1 - exp(-2*pi*j*f/fv))/fv.*sinc(f*Tref);
Hcl_ref = fcw * Hol./(1+Hol);
Hcl_tdc = Hol./(1+Hol);
Hcl_dco = 1./(1+Hol);
Scl_ref = S_ref * abs(Hcl_ref.^2);
Scl_tdc =  PN_tdc .* abs(Hcl_tdc.^2);
Scl_dco = PN_vco .* abs(Hcl_dco.^2);
Scl_tot = Scl_ref + Scl_tdc + Scl_dco;
semilogx(f,10*log10(Scl_ref),'b--');grid on;hold on;
semilogx(f,10*log10(Scl_tdc),'k--');
semilogx(f,10*log10(Scl_dco),'b-');
semilogx(f,10*log10(Scl_tot),'k-');
PN_dcoQ = 1/12.*(KDCO/2^NB_dco./f).^2*1/fref.*sinc(f/fref).^2;
semilogx(f,10*log10(PN_dcoQ.*abs(Hcl_dco.^2)),'b-')

PHE = tckv_div(settle_sample:step_DCO)/Tv*2*pi;
PHE = PHE - mean(PHE);
rbw = 1e4;
[Px_sine,f] = fun_calc_psd_dbs(PHE, fv, rbw);
fstep = f(2)-f(1);
f = f(1:floor(length(f)/2));
Px_sine = Px_sine(1:floor(length(Px_sine)/2));
semilogx(f,Px_sine,'g-');
Px_sine_R = 10.^(Px_sine/10);
sum_Y = sum(Px_sine_R(find(f>1e1)))*fstep;
jitter = sqrt(2*sum_Y) / (2*pi*fv);
xlim([1e3,fv/2]);
ylim([-180,-80]);
title(['fv= ',num2str(fv/1e9),'GHz; rms jitter= ',num2str(floor(jitter*1e15)),' fs']);

%% Spectrum for the ADPLL output
if 0
    Out_sine = sin(2*pi*tckv(1+10000:round(Sim_Time*fcw-2e4)+10000)/Tv);
    figure;plot(Out_sine);
    [Px_sine,f] = fun_calc_psd_dbs(Out_sine, fv, rbw, fstep);
    Px_sine = fftshift(Px_sine);
    Fig = figure('Color',[1 1 1]);
    plot(f+fv/2,Px_sine-10*log10(2));
    grid on;
    xlabel('Frequency (Hz)','fontsize',16);
    ylabel('Phase noise Spectrum (dBc/Hz)','fontsize',16);
    xlim([fv/2 fv*3/2])
    set(gca,'fontsize',16);
    set(findobj(gca,'Type','line'),'LineWidth',2);
    set(gcf,'color','w');
end


%% Present the DTC INL after 2nd-NL cancellation.
figure;
PHRF = 0:1e-3:1;
DCW = ((1-PHRF)*DTC_gain_est )+1;
Delay = DCW*DTC_reso+1*DTC_INL(floor(DCW));
[p,s] = polyfit(PHRF,Delay,1);
Delay_fit = p(1) * PHRF + p(2);
INL_final = Delay - Delay_fit;
plot(PHRF,INL_final,'k-');
hold on;
DCW = ((1-PHRF)*DTC_gain_est + 1*(1-PHRF-0.5).^2*DTC_gain2_est + 1*(1-PHRF-0.5).^3*DTC_gain3_est + 1*(1-PHRF-0.5).^4*DTC_gain4_est + 1*(1-PHRF-0.5).^5*DTC_gain5_est)+1 ;
Delay = DCW*DTC_reso+DTC_INL(floor(DCW));
[p,s] = polyfit(PHRF,Delay,1);
Delay_fit = p(1) * PHRF + p(2);
INL_final = Delay - Delay_fit;
plot(PHRF,INL_final,'r-');