clc;clear;close all;

% FCW_F = 1/16;
% N = 1e5;
% mash_out = 1*order3_mash(FCW_F,N);
% data = mash_out(1e4:end);
% yfft = fft(data)/length(data);
% Yfft = db(abs(yfft));
% fs=125e6;len=length(yfft);
% x = 0:fs/len:fs-fs/len;
% semilogx(x,Yfft,'o-');
% axis tight;
% grid on;
% figure;plot(data,'o-');
% mean(data);

% x = -pi:1e-3:pi;
% y = sin(x);
% plot(x,y,'o-');
% 
% yt = x - x.^3/(3*2)+x.^5/(5*4*3*2) - x.^7/(7*6*5*4*3*2) + x.^9/(9*8*7*6*5*4*3*2) - x.^11/(11*10*9*8*7*6*5*4*3*2) + x.^13/factorial(13);
% 
% hold on;
% plot(x,yt,'*r');
% 
% figure;
% err = yt - y;
% plot(x,err,'o-');