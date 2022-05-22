    function y = order3_mash (x,N)
% x is in the range of [0,1);
% N is the data length of y

dither_en = 1;

mky = zeros(N,1);
y1=1/256;
y2=0;y2_Qprev = 0;
y3=0;y3_Qprev = 0;
y3_diff1 = 0; y3_diff1_prev = 0;

for k = 1:N
    y1 = y1 + x ;%+ dither_en * (randi(2)-1) /2^20;
    y1_Q = floor(y1);
    y1 = y1 - y1_Q;
    
    y2 = y2 + y1 + dither_en * (randi(2)-1) /2^2;
    y2_Q = floor(y2);
    y2 = y2 - y2_Q;
    
    y3 = y3+y2;
    y3_Q = floor(y3);
    y3 = y3-y3_Q;
    
    y2_diff = y2_Q - y2_Qprev;    
    y3_diff1 = y3_Q - y3_Qprev;
    y3_diff2 = y3_diff1 - y3_diff1_prev;

    y2_Qprev = y2_Q;
    y3_Qprev = y3_Q;
    y3_diff1_prev = y3_diff1;

    y(k) = y1_Q + y2_diff + y3_diff2;
    
    Q1(k) = y1_Q;
    Q2(k) = y2_diff;
    Q3(k) = y3_diff2;
 end

end