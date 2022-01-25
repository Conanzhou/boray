%%
k_perp2 = zeros(1,size(yyray,3));
ind = zeros(1,size(yyray,3));
r = zeros(1,size(yyray,3));
z = zeros(1,size(yyray,3));
k_parl = zeros(1,size(yyray,3));
for j=1:size(yyray,3)
    k_perp_temp = yyray(:,4,j).^2 + yyray(:,6,j).^2;
    [k_perp2(j),ind(j)]=min(k_perp_temp);
    r(j) = yyray(ind(j),1,j);
    z(j) = yyray(ind(j),3,j);
    k_parl(j) = yyray(ind(j),5,j) / r(j);
end
k_perp = sqrt(k_perp2);
%%
% figure;
open("./hl2m/output/k_spectrum_2M.fig");
hold on;
[rq,zq] = meshgrid(r,z);
k_rho = diag(interp2(rr',zz',rhorz',rq,zq));

plot(k_rho,k_perp/100,'o','linewidth',2,'MarkerSize',5);
TitleLabels 'k spectrum' 'r(m)' 'k_\perp(cm^{-1})'
% figure;plot3(r,z,k_perp)
% surf(rr,zz,rhorz)
% V = peaks(rq,zq);
% surf(rq,zq,v)