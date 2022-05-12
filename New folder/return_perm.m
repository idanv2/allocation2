function [good_perm,M]=return_perm(v,N,ind2,i0)
%v=0.4;
%N=[0.2;0.2;0.3;0.3];

x=1:1:length(N);
N1=min(N(ind2))-i0(1);
N2=max(N(ind2))-i0(1);
ind_2_x=setdiff(x,ind2);
a=ff2n(size(ind_2_x,2));
M=[];
good_perm=[];
for i=1:1:length(a)
perm=min(N(ind_2_x),v).*(a(i,:)');
if (sum(perm)<=v) && (sum(perm)>=(v-N1-N2))
    good_perm=[good_perm perm];
    m=max(v-sum(perm),0);
    M=[M m];
    %(ind2,ind_2_x,M,good_perm)
end
end
end
%
%for i=1:1:size(M,1)



%end
