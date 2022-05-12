function Q=maxfuncI(x,ind_2,n_2_x,ind_n_2_x,v,M,beta, gamma,epsilon,D,N,k,i0)
n=length(N);
v0=zeros(n,1);
v0(ind_2(1))=x;
v0(ind_2(2))=M-x;
v0(ind_n_2_x)=n_2_x;
%sum(v0);
s0=N-i0-v0;
y0=[s0;v0;i0];
tspan= [0:1e-3:400];opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events', @(t,y) eventfunc(t,y,beta, gamma,epsilon,D,N,k));
[t,y] = ode89(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), tspan, y0,opts);

if t(end)>300
disp('eror')
 [t,y] = ode45(@(t,y) odefcn(t,y,beta, gamma,epsilon,D,N,k), [0:1e-2:300], y0);
I_tot=sum(y(:,2*k+1:3*k),2);
hold on
figure(2)
plot(t,I_tot)
return
end

    
    
   % return

[argvalue, argmax] =max(sum((y(:,2*k+1:3*k)')));
%I_tot=sum(y(:,2*k+1:3*k),2);
hold on
Q=argvalue;
end