clear
load all_dis.dat;
data = all_dis; clear all_dis;

beta = 5;
lambda = 1.8;
[C1, C2] = size(data);
u=data(2:end,:)-lambda*repmat(data(1,:),C1-1,1);
Q=1./(1+exp(beta*u));
Q=sum(Q,2)/C2;
hist(Q,0:0.002:1.0);