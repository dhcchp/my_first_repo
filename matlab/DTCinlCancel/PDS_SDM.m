function [D,R,DMMD] = PDS_SDM (x,N)
% x is in the range of [0,1);
% N is the data length of y
% D is PHRF_D for divider path DTC control
% R is PHRF_R for reference path DTC control

DMMD = zeros(N,1);
R = zeros(N,1);
D = zeros(N,1);

y1=0;
DACC = 0;
URNG1 = unifrnd(0,1,1,N+1);
URNG2 = unifrnd(0,1,1,N);
D2urng = URNG1(2:N+1)-URNG1(1:N);%(URNG1-URNG2);
DQ = 0;
y2_Q = 0;

for k = 1:N
    y1 = y1 + x ;
    y1_Q = floor(y1);
    y1 = y1 - y1_Q;
    
    DACC = DACC + x-y1_Q;
    DQ_prev = DQ;
    y2_Q_prev = y2_Q;
    y2 = DACC + D2urng(k) + DQ_prev; 
    y2_Q = floor(y2);
   % y2 = y2-y2_Q;
    DQ = y2-y2_Q;
    
    DMMD(k) = y1_Q + y2_Q - y2_Q_prev;
    
    DPDS = DACC - y2_Q;
    D(k) = DPDS+D2urng(k);
    R(k) = D2urng(k);

 end

end