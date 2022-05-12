clc
close all
clear
D=[1 1;2 2];
gamma=1;
beta=10;
epsilon=0.1;
k=2;
 v=0.4;
 i0=[0.0001;0.0001];
 v1=v/2;

lb=0;
ub=v;
A=[];
b=[];
Aeq = [];
beq = [];
x = fmincon(@(v1) maxfunc(v1,beta, gamma,epsilon,D),v1,A,b,Aeq,beq,lb,ub)
t=linspace(0,v,100);
Z=t*0;
j=1;
for v1=t
Z(j)=maxfunc(v1,beta, gamma,epsilon,D);
j=j+1;
end
plot(t,Z)
function Q=maxfunc(v1,beta, gamma,epsilon,D)
v=0.3;
v0=[v1;v-v1];
N=[0.6;0.4];
i0=[0.0001;0.0001];
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0 100];opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D), tspan, y0,opts);
%size(t)
%plot(t,y(:,5)+y(:,6))
%hold on
[argvalue, argmax] =max(y(:,5)+y(:,6));
S=[y(argmax,1);y(argmax,2)];
V=[y(argmax,3);y(argmax,4)];
loss=gamma/beta-((D(:,1)')*[S+epsilon*V]);
%scatter(t(argmax),argvalue)
Q=argvalue;
end
function dydt = odefcn(t,y,  beta, gamma,epsilon,D)
S=y(1:2);
V=y(3:4);
I=y(5:6);
ds=-beta*S.*(D*I);
dv=-beta*epsilon*V.*(D*I);
di=beta*(S+epsilon*V).*(D*I)-gamma*I;
dydt=[ds;dv;di];
end
