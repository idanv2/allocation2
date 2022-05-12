function Q=alloc_new(beta,D,N,gamma,epsilon,v)

%load ('GER_data.mat');
%D=contactMatrix;
%N=agDist;
%D=[1 1 1 1 1;2 2 2 2 2;3 3 3 3 3;4 4 4 4 4;5 5 5 5 5];
%N=[0.2;0.2;0.2;0.2;0.2];
%sigma=[1;0.95; 1.975; 2.15;2; 2.05; 2.2; 1.85; 1.85];
%sigma=[1;0.95; 1.975; 2.15; 2.2; 1.85];
%  D=D.*(1./N);
% D=D.*sigma;
%gamma=1;
%epsilon=0.2;
k=length(N);
 %v=0.55;
 i0=zeros(k,1)+1e-8;
 
x=1:1:length(N);
loops=nchoosek(x, 2); % all 2 selections
tot_f=[];
tot_M=[];
tot_v0=[];
tot_loops=[];
tot_solx=[];
tic

parfor i=1:1:length(loops)
ind_2=loops(i,:);
ind_2_x=setdiff((1:1:length(N)),ind_2);
[v0,M]=return_perm(v,N,ind_2,i0);
    for j=1:length(M)
          [solx,f_val]=optim_all(ind_2,ind_2_x,v,M(j),v0(:,j),beta, gamma,epsilon,D,N,k,i0);
          tot_f=[tot_f f_val];
          tot_v0=[tot_v0 v0(:,j)];
          tot_M=[tot_M M(j)];
          tot_solx=[tot_solx solx];
          tot_loops=[tot_loops i];
    end
end   
%save proflie.mat 
toc
%[v_opt]=plot_results(f_min,x_min,perm,N,i0,v);
Q=fig5(tot_f,tot_loops,N,loops,tot_solx,tot_M,tot_v0);




end
