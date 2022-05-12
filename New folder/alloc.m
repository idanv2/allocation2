clc
close all
clear
load ('GER_data.mat');
D=contactMatrix;
N=agDist;
%D=[1 1 1 1 1;2 2 2 2 2;3 3 3 3 3;4 4 4 4 4;5 5 5 5 5];
%N=[0.2;0.2;0.2;0.2;0.2];
%sigma=[1;0.95; 1.975; 2.15;2; 2.05; 2.2; 1.85; 1.85];
sigma=[1;0.95; 1.975; 2.15; 2.2; 1.85];
  D=D.*(1./N);
 D=D.*sigma;
gamma=1;

epsilon=0.2;
k=length(N);
 v=0.55;
 i0=zeros(k,1)+1e-8;

perm=perms(1:k);
beta=5;
min_val=100000;

ind_min=zeros(1,length(perm));
f_min=zeros(1,length(perm));
x_min=zeros(1,length(perm));
tic
parfor i=1:length(perm)
    
    ind=perm(i,:);
   ind_2=(ind(1:2))';
   ind_n_2_x=ind(3:end)';

[solx,f_val]=optim_all(ind_2,ind_n_2_x,v,beta, gamma,epsilon,D,N,k,i0);
f_min(i)=f_val;
x_min(i)=solx;
end
save proflie.mat 
toc
%[v_opt]=plot_results(f_min,x_min,perm,N,i0,v);

function [sol_x,f_val]=optim_all(ind_2,ind_n_2_x,v,beta, gamma,epsilon,D,N,k,i0)
n=length(N);
n_2_x=zeros(1,n-2);
n_2_x(1)=min(N(ind_n_2_x(1))-i0(1),v);
M=max(v-sum(n_2_x),0);
for i=2:n-2
    j=ind_n_2_x(i);
    n_2_x(i)=min(N(j)-i0(1),M);
    M=max(v-sum(n_2_x),0);
end


lb=0;
ub=M;
%ub=N;
A=[];
b=[];
Aeq=[];
beq=[];
x_0=ub/2;
 options=optimset('fmincon');


options=optimset(options,'TolCon',1e-6,'TolX',1e-6,'TolFun',1e-6);
[sol_x,f_val] = fmincon(@(x) maxfuncI(x,ind_2,n_2_x,ind_n_2_x,v,M,beta, gamma,epsilon,D,N,k,i0),x_0,A,b,Aeq,beq,lb,ub,[],options);
end

function Q=maxfuncI(x,ind_2,n_2_x,ind_n_2_x,v,M,beta, gamma,epsilon,D,N,k,i0)
n=length(N);
v0=zeros(n,1);
v0(ind_2(1))=x;
v0(ind_2(2))=M-x;
v0(ind_n_2_x)=n_2_x;
%sum(v0);
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0:1e-3:1000];opts = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @(t,y) eventfunc(t,y,beta, gamma,epsilon,D,N,k));
[t,y] = ode89(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), tspan, y0,opts);

if t(end)>300
disp('eror')
y0
return
end
% [t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), [0:1e-2:300], y0);
%I_tot=sum(y(:,2*k+1:3*k),2);
%hold on
%figure(2)
%plot(t,I_tot)
    
    
   % return

[argvalue, argmax] =max(sum((y(:,2*k+1:3*k)')));
%I_tot=sum(y(:,2*k+1:3*k),2);
hold on
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

