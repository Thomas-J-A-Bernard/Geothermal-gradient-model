clear all; close all;

T0min = -10; T0mid = 0; T0max = 10;
qmmin = 40; qmmid = 50; qmmax = 60; 
q0min = 60; q0mid = 65; q0max = 70;
kmin = 3; k = 3.3; kmax = 3.5;
hr = 10;

y = [0:0.1:30];
F1min = T0min+((qmmin*y)/kmax)+((q0min-qmmin)*hr/kmax)*(1-exp(-y/hr));
F1mid = T0mid+((qmmid*y)/k)+((q0mid-qmmid)*hr/k)*(1-exp(-y/hr));
F1max = T0max+((qmmax*y)/kmin)+((q0max-qmmax)*hr/kmin)*(1-exp(-y/hr));

plot(F1min,y,'m',F1mid,y,'r',F1max,y,'m');
set(gca,'YDir','reverse');



