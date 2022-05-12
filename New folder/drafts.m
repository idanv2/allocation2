M=[1,2,3;4,3,2;1,4,2]
imagesc(M)
for a=1:3
    for b=1:3
        text(b,a,num2str(M(a,b)))
    end
end