% 2021-04-06 15:13 Huasheng XIE, huashengxie@gmail.com, ENN
% read gfile, updated from MA Hao-jie's version
% 21-04-30 23:01 read LiJC's Helicon wave case
% 21-05-02 18:21 modify for b-spline
% 21-05-05 for bi-linear interp2d
close all;clear;clc;
const=1;

if const==1
%     fname='g029300.00800';
%     fname='g36866.860'; 
    fname='g037869.01000';%孔栏位型
%     fname='g037893.01000';%孔栏位型   
else
 %out of region
end   

fid=fopen(fname,'r');
str1=fscanf(fid,'%s',1); str2=fscanf(fid,'%s',1); str3=fscanf(fid,'%s',1); 
str4=fscanf(fid,'%s',1); str5=fscanf(fid,'%s',1);
n3=fscanf(fid,'%d',1);
% n_nomeans=fscanf(fid,'%d',1);%老的hl-2a的nr前有一个0,意义不明
nr=fscanf(fid,'%d',1); nz=fscanf(fid,'%d',1); % nr, nz

Rboxlen=fscanf(fid,'%f',1); Zboxlen=fscanf(fid,'%f',1);
R0=fscanf(fid,'%f',1); Rmin=fscanf(fid,'%f',1); Z0=fscanf(fid,'%f',1);%R0, Z0

Raxis=fscanf(fid,'%f',1); Zaxis=fscanf(fid,'%f',1); Psi_axis=fscanf(fid,'%f',1); %
Psi_bound=fscanf(fid,'%f',1);B0=fscanf(fid,'%f',1);

current=fscanf(fid,'%f',1); Psi_axis1=fscanf(fid,'%f',1); xdum1=fscanf(fid,'%f',1);
Raxis=fscanf(fid,'%f',1); xdum2=fscanf(fid,'%f',1);

Zaxis=fscanf(fid,'%f',1);
xdum3=fscanf(fid,'%f',1);
Psi_bound=fscanf(fid,'%f',1);
xdum4=fscanf(fid,'%f',1);
xdum5=fscanf(fid,'%f',1);

f=fscanf(fid,'%f',nr); %????????????????????

pressure=fscanf(fid,'%f',nr);
ffprime=fscanf(fid,'%f',nr);
pprime=fscanf(fid,'%f',nr);
psi(1:nr,1:nz) =0.0;

% ????gfile????129*129??????????????????????????????????????????????
for i=1:1:nr
    psi(i,1:nz)=fscanf(fid,'%f',nz); %??????????????
end

% %% by zhouyu psi-psi_axis<0,rhorz0 为复数
% psi_min = min(min(psi));
% psi = psi -psi_min;
% Psi_bound = Psi_bound -psi_min;
% Psi_axis = Psi_axis -psi_min;

%%
q=fscanf(fid,'%f',nr); %????????

crop=linspace(Psi_axis,Psi_bound,nr);
R=linspace(Rmin,Rmin+Rboxlen,nr);
% Z=linspace(Z0+Zboxlen*0.5,Z0-Zboxlen*0.5,nz);
Z=linspace(Z0-Zboxlen*0.5,Z0+Zboxlen*0.5,nz);

%????????
dR=R(2)-R(1);
dZ=Z(2)-Z(1);

%??R????????????????????????????????????????????????????????
RR=zeros(nr,nr);
for i=1:nr
    RR(i,1:nr)=R(i);
end

R3=zeros(nr,nz);
for i=1:nr
    R3(1:nz,i)=R(i); %211019,zhouyu
end

%????????Br??Bz
[Bz,Br]=gradient(psi,dR,dZ); %????Br=d(psi)/dz;Bz=d(psi)/dR????????????Br??Bz????????d/dR??d/dZ??????
Br=-Br./R3;
Bz=Bz./R3;
Bp=Br.*Br+Bz.*Bz;
Bp=sqrt(Bp);
[RR2,ZZ2]=meshgrid(R,Z);
% [RR2,ZZ2]=ndgrid(R,Z);

% ????????????
% Psi_axis=-Psi_axis; Psi_bound=-Psi_bound;
% psi_1d=linspace(Psi_axis,Psi_bound,nr);
if Psi_axis>Psi_bound
    psi_1d=linspace(Psi_bound,Psi_axis,nr);
elseif Psi_axis<Psi_bound   
    psi_1d=linspace(Psi_axis,Psi_bound,nr);
else 
    % out of region
end
g=griddedInterpolant(psi_1d,f,'cubic'); %?????????????? x ?????????? v ????????????
gfun=zeros(nr,nz);
Bphi=zeros(nr,nz);
for i=1:nr
    for j=1:nz
        gfun(i,j)=g(psi(i,j));
        Bphi(i,j)=gfun(i,j)./R3(i,j);
    end
end
% %%%%boundary
nbound=fscanf(fid,'%d',1);
nlimiter=fscanf(fid,'%d',1);

RZbound=fscanf(fid,'%f',nbound*2);

RZbound=reshape(RZbound,2,nbound);
Rbound=RZbound(1,:);
Zbound=RZbound(2,:);

RZbound(1,nbound+1)=Raxis;
RZbound(2,nbound+1)=Zaxis;

% limiter%%% 

nbound2=fscanf(fid,'%d',1);
RZlimiter=fscanf(fid,'%f',nlimiter*2);

RZlimiter=reshape(RZlimiter,2,nlimiter);
Rlimiter=RZlimiter(1,:);
Zlimiter=RZlimiter(2,:);

%end limiter %%%%
fclose(fid);


%% 21-04-30 23:06
% q(psi)=dPsi_t/dPsi_p, Psi_p=2*pi*(psi-psi_0)
% Psi_t=0.*psi_1d;

% load('./indexrho=4/helion_ray4.mat');
load('./helion_ray_absorb9.mat');
% load('./helion_hot210504.mat');


if(1==0) % 21-05-04 22:36
    Rtmp=min(R):(max(R)-min(R))/100:max(R);
    Ztmp=min(Z):(max(Z)-min(Z))/200:max(Z);
    [RR2tmp,ZZ2tmp]=meshgrid(Rtmp,Ztmp);
    psitmp=interp2(RR2,ZZ2,psi,RR2tmp,ZZ2tmp,'spline');
    Brtmp=interp2(RR2,ZZ2,Br,RR2tmp,ZZ2tmp,'spline');
    Bztmp=interp2(RR2,ZZ2,Bz,RR2tmp,ZZ2tmp,'spline');
    Bphitmp=interp2(RR2,ZZ2,Bphi,RR2tmp,ZZ2tmp,'spline');
    R=Rtmp; Z=Ztmp; RR2=RR2tmp; ZZ2=ZZ2tmp;
    psi=psitmp; Br=Brtmp; Bz=Bztmp; Bphi=Bphitmp;
    [nr,nz]=size(psi.'); dR=R(2)-R(1); dZ=Z(2)-Z(1);
end

% psilim=psimag+(psilim-psimag)*psifactr;
if(indexrho==2)
    Psi_t=cumsum(q)*(psi_1d(2)-psi_1d(1));
% rho=sqrt(Psi_t/Psi_t(end));
elseif(indexrho==4)
    %     rhorz0=sqrt((psi-Psi_axis)./(Psi_bound-Psi_axis));
    
    Psi_bound1=Psi_bound;
%     Psi_bound1=Psi_axis+(Psi_bound-Psi_axis)*psifactr;
    rhorz0=sqrt((psi-Psi_axis)./(Psi_bound1-Psi_axis));
%     rhorz0=sqrt((psi-psi_min)./(Psi_bound1-psi_min));
end
rhorz=abs(rhorz0); rhorz(rhorz>0.9999)=0.9999;


%% output for ray2d.m ray tracing code
rg=R; zg=Z;
rr=RR2.'; zz=ZZ2.';
% fBr=Br; fBz=Bz;  
fBr=-Br.'; fBz=-Bz.'; % 21-05-01 10:53
fBphi=Bphi.'; fB=sqrt(fBr.^2+fBz.^2+fBphi.^2);
fpsi=psi.';
rhorz=rhorz.';

% rho=rho_bin;
qs0=charge;  qs0(1)=-qs0(1);
ms0=dmas;

n0=1e19;
% genray.dat
if(1==1)
nsrz0(1,:,:)=interp1(rho_bin,squeeze(densprof(:,1)),rhorz)*1e6;
nsrz0(2,:,:)=interp1(rho_bin,squeeze(densprof(:,2)),rhorz)*1e6;
nsrz0(3,:,:)=interp1(rho_bin,squeeze(densprof(:,3)),rhorz)*1e6;
tsrz0(1,:,:)=interp1(rho_bin,squeeze(temprof(:,1)),rhorz)*1e3;
tsrz0(2,:,:)=interp1(rho_bin,squeeze(temprof(:,2)),rhorz)*1e3;
tsrz0(3,:,:)=interp1(rho_bin,squeeze(temprof(:,3)),rhorz)*1e3;
end

nsrz0(isnan(nsrz0))=0.0;
tsrz0(isnan(tsrz0))=0.0;

dr=dR;
dz=dZ;

id0=find(fpsi==min(min(fpsi)));
B0=fB(id0);
R0=rr(id0);
Z0=zz(id0);
% n0=n(id0);

%%
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kgs
qs=qs0*qe;
ms=ms0*me;

fns0=nsrz0;
fts0=tsrz0;


if(1==1)
    
fdBdr=0.*rr; fdBdr(2:(nr-1),:)=(fB(3:nr,:)-fB(1:(nr-2),:))./(2*dr);
fdBdz=0.*rr; fdBdz(:,2:(nz-1))=(fB(:,3:nz)-fB(:,1:(nz-2)))./(2*dz);
fdBrdr=0.*rr; fdBrdr(2:(nr-1),:)=(fBr(3:nr,:)-fBr(1:(nr-2),:))./(2*dr);
fdBrdz=0.*rr; fdBrdz(:,2:(nz-1))=(fBr(:,3:nz)-fBr(:,1:(nz-2)))./(2*dz);
fdBzdr=0.*rr; fdBzdr(2:(nr-1),:)=(fBz(3:nr,:)-fBz(1:(nr-2),:))./(2*dr);
fdBzdz=0.*rr; fdBzdz(:,2:(nz-1))=(fBz(:,3:nz)-fBz(:,1:(nz-2)))./(2*dz);
fdBphidr=0.*rr; fdBphidr(2:(nr-1),:)=(fBphi(3:nr,:)-fBphi(1:(nr-2),:))./(2*dr);
fdBphidz=0.*rr; fdBphidz(:,2:(nz-1))=(fBphi(:,3:nz)-fBphi(:,1:(nz-2)))./(2*dz);

fdns0dr=0.*fns0; fdns0dr(:,2:(nr-1),:)=(fns0(:,3:nr,:)-fns0(:,1:(nr-2),:))./(2*dr);
fdns0dz=0.*fns0; fdns0dz(:,:,2:(nz-1))=(fns0(:,:,3:nz)-fns0(:,:,1:(nz-2)))./(2*dz);
else  % 21-04-07 13:27 fixed r,z index, which seems should be opposite
fdBrdr=0.*rr; fdBrdr(:,2:(nr-1))=(fBr(:,3:nr)-fBr(:,1:(nr-2)))./(2*dr);
fdBrdz=0.*rr; fdBrdz(2:(nz-1),:)=(fBr(3:nz,:)-fBr(1:(nz-2),:))./(2*dz);
fdBzdr=0.*rr; fdBzdr(:,2:(nr-1))=(fBz(:,3:nr)-fBz(:,1:(nr-2)))./(2*dr);
fdBzdz=0.*rr; fdBzdz(2:(nz-1),:)=(fBz(3:nz,:)-fBz(1:(nz-2),:))./(2*dz);
fdBphidr=0.*rr; fdBphidr(:,2:(nr-1))=(fBphi(:,3:nr)-fBphi(:,1:(nr-2)))./(2*dr);
fdBphidz=0.*rr; fdBphidz(2:(nz-1),:)=(fBphi(3:nz,:)-fBphi(1:(nz-2),:))./(2*dz);
fdBdr=0.*rr; fdBdr(:,2:(nr-1))=(fB(:,3:nr)-fB(:,1:(nr-2)))./(2*dr);
fdBdz=0.*rr; fdBdz(2:(nz-1),:)=(fB(3:nz,:)-fB(1:(nz-2),:))./(2*dz);

fdnedr=0.*rr; fdnedr(:,2:(nr-1))=(fne(:,3:nr)-fne(:,1:(nr-2)))./(2*dr);
fdnedz=0.*rr; fdnedz(2:(nz-1),:)=(fne(3:nz,:)-fne(1:(nz-2),:))./(2*dz);
fdnidr=0.*rr; fdnidr(:,2:(nr-1))=(fni(:,3:nr)-fni(:,1:(nr-2)))./(2*dr);
fdnidz=0.*rr; fdnidz(2:(nz-1),:)=(fni(3:nz,:)-fni(1:(nz-2),:))./(2*dz);
% fdnedr=0.*rr; fdnedr(2:(nr-1),:)=(fne(3:nr,:)-fne(1:(nr-2),:))./(2*dr);
% fdnedz=0.*rr; fdnedz(:,2:(nz-1))=(fne(:,3:nz)-fne(:,1:(nz-2)))./(2*dz);
% fdnidr=0.*rr; fdnidr(2:(nr-1),:)=(fni(3:nr,:)-fni(1:(nr-2),:))./(2*dr);
% fdnidz=0.*rr; fdnidz(:,2:(nz-1))=(fni(:,3:nz)-fni(:,1:(nz-2)))./(2*dz);
end

S=length(qs);
save('genray_eqdata.mat','rg','zg','dr','dz','rr','zz',...
    'fB','fBr','fBz','fBphi','fns0','fts0','fdBdr','fdBdz',...
    'fdBrdr','fdBrdz','fdBzdr','fdBzdz','fdBphidr','fdBphidz',...
    'fdns0dr','fdns0dz','qs','ms','R0','Z0','B0','n0','fpsi','S');

%% plot
h=figure('unit','normalized','Position',[0.01 0.05 0.8 0.6],...
  'DefaultAxesFontSize',14);
subplot(241);
contour(rr,zz,fpsi,100);
xlabel('R');ylabel('Z');title('\psi');
axis equal; colorbar;

subplot(242);
contour(rr,zz,fB,100);
xlabel('R');ylabel('Z');title('B');
axis equal; colorbar;

subplot(243);
contour(rr,zz,fBr,100);
xlabel('R');ylabel('Z');title('B_r');
axis equal; colorbar;

subplot(244);
contour(rr,zz,fBz,100);
xlabel('R');ylabel('Z');title('B_z');
axis equal; colorbar;

subplot(245);
contour(rr,zz,fBphi,100);
xlabel('R');ylabel('Z');title('B_t');
axis equal; colorbar;

subplot(246);
contour(rr,zz,squeeze(fns0(1,:,:)),100);
xlabel('R');ylabel('Z');title('n_{e}');
axis equal; colorbar;

subplot(247);
contour(rr,zz,squeeze(fts0(1,:,:)),100);
xlabel('R');ylabel('Z');title('T_{e}');
axis equal; colorbar;

subplot(248);
contour(rr,zz,squeeze(fdns0dz(1,:,:)),100);
xlabel('R');ylabel('Z');title('dn_{s1}/dz');
axis equal; colorbar;

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['gfile_eqdata_B0=',num2str(B0,3),...
    'T,ni0=',num2str(n0,3),'.png']);
%%

f=freqcy;
w=2*pi*f;
c=2.99792458e8;
wkr=wn_r/c*w;
wkz=wn_z/c*w;
wkphi=wn_phi/c*w;
br=sb_r*1e-4; bz=sb_z*1e-4; bphi=sb_phi*1e-4; btot=sbtot*1e-4;
ne=sene*1e6; te=ste*1e3;
wr=wr/100; wz=wz/100;
vgr_r=vgr_r*c; vgr_phi=vgr_phi*c; vgr_z=vgr_z*c;
flux_r=flux_r*c; flux_phi=flux_phi*c; flux_z=flux_z*c; 
delpwr=delpwr*1e-7;
% salphas=
squeeze(salphas(:,id,1));

id=1;
yyid=[wr(:,id),wphi(:,id),wz(:,id),wkr(:,id),wkphi(:,id).*wr(:,id),wkz(:,id),...
    wnpar(:,id)*w/c,wnper(:,id)*w/c,br(:,id),bz(:,id),bphi(:,id),btot(:,id),...
    ne(:,id),te(:,id),vgr_r(:,id),vgr_phi(:,id)./wr(:,id),vgr_z(:,id),...
    flux_r(:,id),flux_phi(:,id)./wr(:,id),flux_z(:,id),delpwr(:,id)];
save(['rayid=',num2str(id),'.mat'],'yyid','f','iabsorp');

% id=1;

yyallray=zeros(size(yyid,1),size(yyid,2),size(wr,2));
for id=1:size(wr,2)
    yyallray(:,:,id)=yyid;
    save('ray_all.mat','yyallray','f','iabsorp');
end

%% check the density and temperature profiles
if(1==1)
figure;
rtmp=w_r_densprof_nc; nstmp=w_dens_vs_r_nc*1e13; Tstmp=w_temp_vs_r_nc*1e3;

izc=floor(nz/2)+1;
subplot(121);
plot(rtmp,nstmp(:,1),'linewidth',2); hold on;
plot(rr(:,izc),squeeze(fns0(1,:,izc)),'r--','linewidth',2);
xlabel('R'); ylabel('n_{s0}');

subplot(122);
plot(rtmp,Tstmp(:,1),rr(:,izc),squeeze(fts0(1,:,izc)),'r--','linewidth',2);
xlabel('R');ylabel('T_s');
end
