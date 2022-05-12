function [value,isterminal,direction]=eventfunc(t,y,beta, gamma,epsilon,D,N,k)
S=y(1:k);
V=y(k+1:2*k);
I=y(2*k+1:3*k);
di=beta*(S+epsilon*V).*(D*(I./N))-gamma*I;
value=sum(di);
isterminal=1;
direction=-1;
end
