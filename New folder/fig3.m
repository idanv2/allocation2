clc
close all
clear

D=[1 1;2 2];
N=[0.6;0.4];
beta=0.3;
gamma=1;
epsilon=0.1;
k=length(N);
 v=0.1;
 v1=0.05;
 i0=zeros(k,1)+1e-6;
%t=linspace(0,1,100);
%f=root2d(t,beta, gamma,epsilon,D,N,k,i0,v1,v);
%plot(t,f)
%criteria
bt=beta/gamma;
v2=v-v1;
S1=N(1)-i0(1)-v1;
S2=N(2)-i0(2)-v2;
e1=(D(1,1)*(S1+epsilon*v1)+D(2,2)*(S2+epsilon*v2));
e2=D(1,1)*(S1*exp(-bt*D(1,1))+epsilon*v1*exp(-bt*epsilon*D(1,1)))+...
    D(2,2)*(S2*exp(-bt*D(2,2))+epsilon*v2*exp(-bt*epsilon*D(2,2)));
clc
gamma/beta-e1
gamma/beta-e2
%%%
v1=linspace(0,v,20);
H=[];
for i=1:20
z=sumR(beta, gamma,epsilon,D,N,k,i0,v1(i),v);
H=[H matara(z,beta, gamma,epsilon,D,N,k,i0,v1(i),v)];
end
plot(v1,H)
function Z=sumR(beta, gamma,epsilon,D,N,k,i0,v1,v)
Z = fsolve(@(z) root2d(z,beta, gamma,epsilon,D,N,k,i0,v1,v),0);

end

function H = root2d(z,beta, gamma,epsilon,D,N,k,i0,v1,v)
bt=beta/gamma;
v2=v-v1;
x=exp(-bt*z);
H=D(1,1)*(N(1)-v1-i0(1))*x.^(D(1,1))+...
    D(1,1)*epsilon*v1*x.^(epsilon*D(1,1))+...
    D(2,2)*(N(2)-v2-i0(2))*x.^(D(2,2))+...
    D(2,2)*epsilon*v2*x.^(epsilon*D(2,2));
H=H-gamma/beta;

end
function H = matara(z,beta, gamma,epsilon,D,N,k,i0,v1,v)
bt=beta/gamma;
v2=v-v1;
x=exp(-bt*z);
H=z+(N(1)-v1-i0(1))*x.^(D(1,1))+...
    v1*x.^(epsilon*D(1,1))+...
    (N(2)-v2-i0(2))*x.^(D(2,2))+...
    v2*x.^(epsilon*D(2,2));


end










