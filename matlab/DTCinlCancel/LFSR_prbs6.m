clc;clear; close all;

N = 1e5;
Q = zeros(1,6);
y = zeros(1,N);
Q(1) = 0;
for k = 1:N
    D1 = 1-xor(Q(5),Q(6));
%     D1 = xor(D1,Q(5));
%     D1 = xor(D1,Q(4));
    Q = [D1,Q(1:5)];
    y(k) = Q(6);
end

yfft = fft(y);
Yfft = db(abs(yfft));

semilogx(Yfft,'o-');
grid on;

