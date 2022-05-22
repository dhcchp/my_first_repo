clc;clear;close all;

FCW_F = 0.1111+1/16;
N = 1e5;
mash_out = 1*order3_mash(FCW_F,N);
data = mash_out(1:end);
yfft = fft(data)/length(data);
Yfft = db(abs(yfft));
fs=125e6;len=length(yfft);
x = 0:fs/len:fs-fs/len;
semilogx(x,Yfft,'o-');
axis tight;
grid on;
figure;plot(data,'o-');
mean(data);