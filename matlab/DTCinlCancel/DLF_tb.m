clc;clear;close all;
P=bodeoptions;
P.FreqUnits='Hz';

s = tf('s');
H_filter = 2^(-1) + 2^(-3)/s;
bode(H_filter,P);
grid on;

Hol = H_filter * 1/s;
Hcl = Hol/(1+Hol);
figure;
bode(Hcl,P);
grid on;

len = 1e5;
x = rand(len,1);
filter_in = zeros(size(x));
filter_in_accum = zeros(size(x));
filter_out = zeros(size(x));
vco_out = zeros(size(x));

num = H_filter.num{1};
alpha = num(1);rho = num(2);

for k = 2:len
    filter_in(k) = x(k) - vco_out(k-1);
    filter_in_accum (k) = filter_in_accum(k-1) + filter_in(k);
    filter_out(k) = alpha * filter_in(k) + rho * filter_in_accum(k);
    vco_out(k) = vco_out(k-1) + filter_out(k);
end

yfft = fft(vco_out);
Yfft = db(abs(yfft));
freq = 0:2/len:2-2/len;
semilogx(freq,Yfft,'o-');grid on;