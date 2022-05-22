function y = order2_mash (x,N)
% x is in the range of [0,1);
% N is the data length of y

dither_en = 0;

y = zeros(N,1);
y1=0;
y2=0;y2_Qprev = 0;

for k = 1:N
    y1 = y1 + x + dither_en * (randi(2)-1) /2^10;
    y1_Q = floor(y1);
    y1 = y1 - y1_Q;
    
    y2 = y2 + y1;
    y2_Q = floor(y2);
    y2 = y2 - y2_Q;

    y2_diff = y2_Q - y2_Qprev;    

    y2_Qprev = y2_Q;

    y(k) = y1_Q + y2_diff;
 end

end