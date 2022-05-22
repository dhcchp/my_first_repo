%% General Description
% This is modified from the model for J Choi ISSCC 2021, TIPM architecture
% Note that DCWR stands different compared with ISSCC 2020 version. 
% A 12'bit DTC is used. Consider it as 12 separate DTCs, calibrate the
% binary bits independently.

% Version 1.0, Peng CHEN, 2021.12.08
% Version 1.1, Peng CHEN, 2022.05.17, This version fails. Maybe NerualNetwork can
% help? The appearance of "12 DTCs" may need reshape.

%% Initial Parameter Settings
clc;clear;close all;
Sim_Time = 2^13; % The number is in FREF cycles.
settle_sample = 6e6; % The unit is in CKV cycles.
rng(23);

%% DTC alone test
FCW_F = 0.11;
B = 12;
dcw_r = zeros(1,B);
DCW_R = zeros(1,B);
track_g = zeros(Sim_Time,B);
PHRF = FCW_F * [1:Sim_Time]+1/2^10*rand(1,Sim_Time);
PHRF = 1*rand(1,Sim_Time);
OV = zeros(1,Sim_Time);
tmp = zeros(size(OV));
PHRF = mod(PHRF,1);
%plot(PHRF,'o-');

% w = [1./2.^[1:B],1/2^B];
% g = ones(1,B+1);
w = [1./2.^[1:B]];
w_t = 1e-12*[2048,1024,512,256,128,64,32,16,8,4,2,1] * 1e3/4096;
g = ones(1,B);
hist_DCWR = zeros(1,B);
mu = -1e-1;
KDTC=0.8;
for k = 1:Sim_Time
    dcw_r = zeros(1,B);
    DCW_R = zeros(1,B);
    dcw_s = dec2bin(PHRF(k)*2^B*0.7);
    for k2 = 1:length(dcw_s)
        dcw_r(k2+B-length(dcw_s)) = str2double(dcw_s(k2));
    end
    DCW_r = sum(dcw_r.*w.*g);
    if DCW_r >= 1-1/2^B
        DCW_r = 1-1/2^B;
    end
    if DCW_r < 0
        DCW_r = 0;
    end
    DCW_s = dec2bin(DCW_r*2^B);
    for k2 = 1:length(DCW_s)
        DCW_R(k2+B-length(DCW_s)) = str2double(DCW_s(k2));
    end
    tdelay = sum(DCW_R.*w_t);
%     phe = floor((tdelay - PHRF(k)*1000e-12*KDTC)/1e-12);
    phe = ((tdelay - PHRF(k)*1000e-12*0.7)/1e-12);
    for k3 = 1:B
%         g(k3) = g(k3) + mu * phe * 1/2^(B-k3) * (DCW_R(k3)-0.5);
%         track_g(k,k3) = g(k3);
          track_g(k,k3) =  mu * phe * 1/2^(B-k3) * (DCW_R(k3)-0.5);
    end
%     KDTC = KDTC-mu*phe;
    tmp(k) = phe;
    hist_DCWR = hist_DCWR+DCW_R;
end

plot(track_g(:,1:12),'o-');
figure;
plot(tmp,'o-');
figure;plot(hist_DCWR,'o-');
sum(track_g)*100