% load('yyray.mat', 'yyray');
% figure();hold on;
%%
num_x = length(yyray(:,1,1));
[point_x,point_o] = deal(ones(nray,1));
for j=2:nray
    [L,ind2]=deal(ones(num_x,1));
    for i=1000:num_x
        [L(i),ind2(i)]=leng([yyray(i,1,1),yyray(i,3,1)],[yyray(:,1,j),yyray(:,3,j)]);
    end
    [~,point_x(j)] = min(L(:,1));
    point_o(j) = ind2(point_x(j));
    quiver(yyray(point_x(j),1,1),yyray(point_x(j),3,1),yyray(point_x(j),4,1)*2e-5,yyray(point_x(j),6,1)*2e-5,'b','LineWidth',2);
    quiver(yyray(point_o(j),1,j),yyray(point_o(j),3,j),yyray(point_o(j),4,j)*2e-5,yyray(point_o(j),6,j)*2e-5,'r','LineWidth',2);
end

%%
function [len,ind]=leng(p,l)
    L=sum( (l-p).^2 ,2);
    [len,ind]=min(L);
end