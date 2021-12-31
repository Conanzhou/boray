% x1x,x1y,o1x,o1y
[L,ind2]=deal(ones(length(x1x),1));
for i=1000:length(x1x)
    [L(i),ind2(i)]=leng([x1x(i),x1y(i)],[o1x,o1y]);
end
[~,ind]=min(L(:,1));

%%
L1=x1x+x1y*1i;
L2=o1x+o1y*1i;

Length=abs(L1(1000:end)'-L2(1500:end));
[~,ind]=min(Length(:));

%%
function [len,ind]=leng(p,l)
L=sum( (l-p).^2 ,2);
[len,ind]=min(L);
end

