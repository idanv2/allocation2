clc
close all
clear
load ('USA_data.mat');
D=contactMatrix;
N=agDist;

%D=[1 2 3;3 4 5;2 1 1];
D=[1 1;2 2];
N=[0.6;0.4];
%D=D/max(eigs(D));
%N=[0.3;0.3;0.4];
%D=[1 1 1;2 2 2;3 3 3];
%N=[0.3;0.4;0.3];

%D=D/max(eigs(D));
gamma=1;
epsilon=0.1;
k=length(N);
 v=0.1;
 i0=zeros(k,1)+1e-6;
 v0=zeros(k,1)+v/k;
lb=zeros(k,1);
ub=N-i0;
%ub=N;
A=[];
b=[];
Aeq = ones(1,k);
beq = v;
B=[];
C=[];

for beta=[6:0.25:10]
x = fmincon(@(v0) maxfuncI(v0,beta, gamma,epsilon,D,N,k,i0),v0,A,b,Aeq,beq,lb,ub);
%y=fmincon(@(v0) maxfuncI3(v0,beta, gamma,epsilon,D,N,k,i0),v0,A,b,Aeq,beq,lb,ub);
B=[B x/v];
%C=[C y/v];
end
v0=x;
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0 50];opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), tspan, y0,opts);
I_tot=sum(y(:,2*k+1:3*k),2);
plot (t,I_tot)

plotgraph=1;

if plotgraph==1
subplot(1,2,1)
imagesc(B);
colormap;
colorbar;
yticks([1:k]);
%xticks([2:20]);
xlabel('\beta');
ylabel('$$\vec{v}_{optimal}$$','interpreter','latex');
 hYLabel = get(gca,'YLabel');
 set(hYLabel,'rotation',0,'VerticalAlignment','middle');
title('Optimal allocation','interpreter','latex');

axis square

subplot(1,2,2)
imagesc(D);
set(gca,'xtick',[])
set(gca,'ytick',[])
title('Contact matrix','interpreter','latex')
yticks([1:k])
xticks([1:k])
colorbar
axis square


end
function Q=maxfuncI3(v0,beta, gamma,epsilon,D,N,k,i0)
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0 200];opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), tspan, y0,opts);
%[argvalue, argmax] =max(sum((y(:,3*k-2:3*k)')));
[argvalue, argmax] =max((y(:,3*k)'));
Q=argvalue;
end
function Q=maxfuncI(v0,beta, gamma,epsilon,D,N,k,i0)
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0 100];opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), tspan, y0,opts);
%size(t)
%plot(t,sum((y(:,2*k+1:3*k)')))
%hold on
[argvalue, argmax] =max(sum((y(:,2*k+1:3*k)')));
S=y(argmax,1:k)';
V=y(argmax,k+1:2*k)';
I=y(argmax,2*k+1:3*k)';
I_tot=sum(y(:,2*k+1:3*k),2);
if sum(y(end,2*k+1:3*k),2)>1e-6
    close all
     plot(t,I_tot)
    clc
disp('error')
return 
end
%figure(1)
%hold on
%plot(I_tot)
%loss=sum(beta*(S+epsilon*V).*(D*I)-gamma*I);
%loss
%loss=gamma/beta-((D(:,1)')*(S+epsilon*V));
%scatter(t(argmax),argvalue)
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

