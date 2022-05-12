function [sol_x,f_val]=optim_all(ind_2,ind_n_2_x,v,M,v0,beta, gamma,epsilon,D,N,k,i0)
lb=max(0,M-N(ind_2(2))+i0(1) );
ub=min(M,N(ind_2(1))-i0(1) ) ;

%A=[];
%b=[];
%Aeq=[];
%beq=[];
%x_0=ub/2;
 %options=optimset('fmincon');
%options=optimset(options,'TolCon',1e-6,'TolX',1e-6,'TolFun',1e-6);
%[sol_x,f_val] = fmincon(@(x) maxfuncI(x,ind_2,v0,ind_n_2_x,v,M,beta, gamma,epsilon,D,N,k,i0),x_0,A,b,Aeq,beq,lb,ub,[],options);
[sol_x,f_val]=fminbnd(@(x) maxfuncI(x,ind_2,v0,ind_n_2_x,v,M,beta, gamma,epsilon,D,N,k,i0),lb,ub);
end