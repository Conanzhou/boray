        fig_rz=figure('unit','normalized','Position',[0.01 0.05 0.3 0.6],'DefaultAxesFontSize',14);hold on;
        ylim([-1.2 1.2]);
%         contour(rr,zz,fpsi,100);
%         plot(yy(1,1),yy(1,3),'rx','linewidth',2);hold on;
%         contour(rr,zz,rhorz,'LevelList',0:0.05:0.7);hold on;
%         contour(rr,zz,rhorz,'LevelList',0.7:0.01:1.1);
%         contour(rr,zz,squeeze(fns0(1,:,:)),100);
        contour(rr,zz,squeeze(fns0(1,:,:)),'LevelList',(0.6:0.01:2)*1e19);
        xlabel('R'); ylabel('Z');box on;
        
        num_theta = 1;
        num_f = 7;
        for i = 1:num_f
            for j = 1:num_theta
                temp_num = (i-1)*num_theta+j;
                plot(yyray(:,1,temp_num),yyray(:,3,temp_num),'.');
            end
        end