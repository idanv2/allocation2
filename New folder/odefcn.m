function dydt = odefcn(t,y,  beta, gamma,epsilon,D,N,k)

S=y(1:k);
V=y(k+1:2*k);
I=y(2*k+1:3*k);
ds=-beta*S.*(D*(I./N));
dv=-beta*epsilon*V.*(D*(I./N));
di=beta*(S+epsilon*V).*(D*(I./N))-gamma*I;
dydt=[ds;dv;di];
end