function Q=fig5(tot_f,tot_loops,N,loops,tot_solx,tot_M,tot_v0)
j=find(tot_f==min(tot_f));
opt_2=loops(tot_loops(j),:);
opt_rest=setdiff((1:1:length(N)),opt_2);
v0=zeros(length(N),1);
v0(opt_2)=[tot_solx(j);tot_M(j)-tot_solx(j)];
v0(opt_rest)=tot_v0(:,j);
Q=v0';
end


