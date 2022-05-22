function y = order1_mash (x,N)
% x is in the range of [0,1);
% N is the data length of y

dither_en = 1;

y = zeros(N,1);
y1=0;

for k = 1:N
    y1 = y1 + x + dither_en * (randi(2)-1) /2^4;
    y1_Q = floor(y1);
    y1 = y1 - y1_Q;
    
    y(k) = y1_Q;
 end

end