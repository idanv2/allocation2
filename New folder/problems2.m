clc
close all
clear
D=[1 1;2 2];
N=[0.6;0.4];
D=D/max(eigs(D));
gamma=1;
epsilon=0.1;
k=length(N);
 v=0.2;
 
 i0=zeros(k,1)+1e-5;
 v0=zeros(k,1)+v/k;
lb=zeros(k,1);
ub=N-i0;
%ub=N;
A=[];
b=[];
Aeq = ones(1,k);
beq = v;
B=[];
C=D(:,1);
for beta=[17]
x = fmincon(@(v) maxfuncI(v,beta, gamma,epsilon,C,N,k,i0),v0,A,b,Aeq,beq,lb,ub);
B=[B x/v];
end
%x = fmincon(@(v0) maxfuncI(v0,beta, gamma,epsilon,C,N,k,i0),v0,A,b,Aeq,beq,lb,ub);
function F=maxfuncI(v0,beta, gamma,epsilon,C,N,k,i0)
z0=0.1;
F=fsolve(@(z) Z(z,v0,beta, gamma,epsilon,C,N,k,i0),z0);
%Z(F,v0,beta, gamma,epsilon,C,N,k,i0)
tspan= [0 50];opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
s0=N-i0-v0;
y0=[s0;v0;i0;i0*0];
[t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,[C C],N,k), tspan, y0,opts);
[argvalue, argmax] =max(sum((y(:,2*k+1:3*k)')));
R_tot=sum(y(argmax,3*k+1:4*k)');
%plot(t,R_tot)
end
function x=Z(z,v0,beta, gamma,epsilon,C,N,k,i0)
s0=N-i0-v0;
E1=s0.*exp(-beta*C*z);
E2=epsilon*v0.*exp(-beta*epsilon*C*z);
x=sum(C.*(E1+E2))-gamma/beta;
end
function dydt = odefcn(t,y,  beta, gamma,epsilon,D,N,k)

S=y(1:k);
V=y(k+1:2*k);
I=y(2*k+1:3*k);
ds=-beta*S.*(D*I);
dv=-beta*epsilon*V.*(D*I);
di=beta*(S+epsilon*V).*(D*I)-gamma*I;
dr=gamma*I;
dydt=[ds;dv;di;dr];
end
