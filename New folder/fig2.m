clc
close all
clear
load ('USA_data.mat');
D=contactMatrix*100;
N=agDist;
D(6:end,:)=[];
D(:,6:end)=[];
N(6:end)=[];
gamma=1;
epsilon=0.1;
k=length(N);
 v=0.2;
 i0=zeros(k,1)+1e-6;
 v0=zeros(k,1)+v/k; % initial guess
v0(3)=0;
v0(4)=0;
v0(2)=0;
v0(1)=0;
v0(5)=0.2;


beta=2;
a1=0.1;
a2=0.05;
clc
v0=[  a1    a2   0.0000    0.0000    0.2-a1-a2]'
maxfuncI(v0,beta, gamma,epsilon,D,N,k,i0)
function Q=maxfuncI(v0,beta, gamma,epsilon,D,N,k,i0)
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0:1e-4:5];opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @(t,y) eventfunc(t,y,beta, gamma,epsilon,D,N,k));
[t,y] = ode89(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), tspan, y0,opts);
%size(t)
%plot(t,sum((y(:,2*k+1:3*k)')))
%hold on
[argvalue, argmax] =max(sum((y(:,2*k+1:3*k)')));
[t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), [0:1e-2:10], y0);
I_tot=sum(y(:,2*k+1:3*k),2);

Q=argvalue;
end
function dydt = odefcn(t,y,  beta, gamma,epsilon,D,N,k)

S=y(1:k);
V=y(k+1:2*k);
I=y(2*k+1:3*k);
ds=-beta*S.*(D*I);
dv=-beta*epsilon*V.*(D*I);
di=beta*(S+epsilon*V).*(D*I)-gamma*I;
dydt=[ds;dv;di];
end
function [value,isterminal,direction]=eventfunc(t,y,beta, gamma,epsilon,D,N,k)


S=y(1:k);
V=y(k+1:2*k);
I=y(2*k+1:3*k);
di=beta*(S+epsilon*V).*(D*I)-gamma*I;
value=sum(di);
isterminal=1;
direction=-1;
end










