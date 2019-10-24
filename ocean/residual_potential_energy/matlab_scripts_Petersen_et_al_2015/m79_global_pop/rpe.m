% script rpe
% compute resting potential energy, ala Ilicak ea 2012 for POP data
% Mark Petersen, LANL, Dec 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%load '~/mpas/email/130911_Ilicak_fig_data/fig17.mat'

cl=[...
  0    0.6  0.9 ; ... % blue
  1    0    0   ; ... % red
  0.1  1    0.1 ; ... % green
  0    0    0   ; ... % black
  0    0.4  0.1 ; ... % dark green
  1    0.4  0   ; ... % orange
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];
marker='*s^v';
title_txt='bc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10 km 100 day


dir={'m79'}; abc='bcda'; dims=1*[16 56];

%nu_h=[1 5 10 20 200];
nu_h=[20];
grid_spacing=10e3;

hist_files = {...
't.x1_RPE2a.00000101-00000312all.nc'
't.x1_RPE3a.00000101-00000312all.nc'
't.x1_RPE4a.00000101-00000312all.nc'
't.x1_RPE1a.00000101-00000312all.nc'
    };

time_fields=[1:14];
min_n=1;max_n = length(time_fields);

kmt_file = 'topography_20010702.ieeei4';
grid_file = 'horiz_grid_20010402.ieeer8';
depth_file = 'in_depths.dat';

title_txt={
    'POP 1 degree',...
	  };
if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt10km100day = zeros(1,length(abc));
keMeanTime10km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  DrpeDtPOPx1 = zeros(length(time_fields),length(abc));
nVertLevels = 40;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
[time,rpeTot(:,j),rpePert(:,j),rpeNorm(:,j),DrpeDtPOPx1(:,j),meanDrpeDt10km(j),keMeanTime10km(j)] ...
    = sub_rpe(wd,dir,abc(j),hist_files(j), kmt_file, grid_file, depth_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt10km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt10km100day(j));
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
%  keMeanTime(j) = mean(keMeanVolume(min_n:max_n));
%  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
%  DrpeDtNorm = derivative_ord2(time,rpeNorm);
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
end

DrpeDtNorm = zeros(size(rpeNorm));
for j=1:length(abc)
%  DrpeDtNorm(:,j) = derivative_ord2(time,squeeze(rpeNorm(:,j))');
end

end

meanDrpePertDtAll = zeros(4,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dRPE/dt with time POP 1deg
load('data/popx1_140127.mat')

figure(14);clf 
nsims=length(abc);
for j=1:nsims
%  h=semilogy((time-time(1))/3600/24,DrpePertDt(:,j),'-.');
  h=semilogy((time-time(1))/3600/24,DrpePertDt(:,j),'-.');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
  meanDrpePertDtAll(3,j) = mean(DrpePertDt(5:length(time),j));
end

size(DrpePertDt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dRPE/dt with time MPAS 120km

load('data/m55edcb_global_spindown.mat')
DrpePertDt = zeros(size(rpePert));
for j=1:length(abc)
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
end

figure(14);
nsims=length(abc);
for j=1:nsims
  h=semilogy((time-time(1))/3600/24,DrpePertDt(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
  meanDrpePertDtAll(1,j) = mean(DrpePertDt(21:length(time),j));
end

size(DrpePertDt)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dRPE/dt with time MPAS 60km

load('data/m55jihk_global_spindown.mat')
DrpePertDt = zeros(size(rpePert));
for j=1:length(abc)
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
end
figure(14);
nsims=length(abc);
for j=1:nsims
  h=semilogy((time-time(1))/3600/24,DrpePertDt(:,j),'--');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
  meanDrpePertDtAll(2,j) = mean(DrpePertDt(21:length(time),j));
end
size(DrpePertDt)

%%%%%%%%%%%% add Ilicak data
load '~/mpas/email/130911_Ilicak_fig_data/fig20.mat'

j=1;
h=semilogy(mom_kd_1e_5_time,mom_kd_1e_5_drpe_dt,':');
  meanDrpePertDtAll(4,1) = mean(mom_kd_1e_5_drpe_dt(20:69));
set(h,'Color',cl(j,:),'LineWidth',1); j=j+1;

h=semilogy(mom_kd_1e_6_time,mom_kd_1e_6_drpe_dt,':');
  meanDrpePertDtAll(4,2) = mean(mom_kd_1e_6_drpe_dt(20:69));
set(h,'Color',cl(j,:),'LineWidth',1); j=j+1;

h=semilogy(mom_kd_1e_7_time,mom_kd_1e_7_drpe_dt,':');
  meanDrpePertDtAll(4,3) = mean(mom_kd_1e_7_drpe_dt(20:69));
set(h,'Color',cl(j,:),'LineWidth',1); j=j+1;

h=semilogy(mom_kd_0_time,mom_kd_0_drpe_dt,':');
  meanDrpePertDtAll(4,4) = mean(mom_kd_0_drpe_dt(20:69));
set(h,'Color',cl(j,:),'LineWidth',1); 


grid on
xlabel('time, days')
ylabel('dRPE/dt, W/m^2')
axis([0 70 8e-4 1e-2])
group_legend={...
'POP 1^o \kappa_v=1e-5','POP 1^o \kappa_v=1e-6','POP 1^o \kappa_v=1e-7','POP 1^o \kappa_v=0',...
'MPAS-O 120km \kappa_v=1e-5','MPAS-O 120km \kappa_v=1e-6','MPAS-O 120km \kappa_v=1e-7','MPAS-O 120km \kappa_v=0',...
'MPAS-O 60km \kappa_v=1e-5','MPAS-O 60km \kappa_v=1e-6','MPAS-O 60km \kappa_v=1e-7','MPAS-O 60km \kappa_v=0',...
'MOM 1^o \kappa_v=1e-5','MOM 1^o \kappa_v=1e-6','MOM 1^o \kappa_v=1e-7','MOM 1^o \kappa_v=0',...
	      };

%group_legend={...
%'120km \nu_h=1e12','120km \nu_h=4e12','120km \nu_h=1e13','120km \nu_h=4e13','120km \nu_h=1e14',...
%'60km \nu_h=4e11','60km \nu_h=1e12','60km \nu_h=4e12','60km \nu_h=1e13','60km \nu_h=4e13',...
%	      };
legend(group_legend,'location','EastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 5 ])
fig=[char(dir) abc(1:nsims) '_DrpePertDtPOP'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);


%%%%%%%%%%%% bar chart

figure(21); clf;
h=bar(meanDrpePertDtAll,'ShowBaseLine', 'off')
set(h,'BaseValue',1e-4);
set(gca,'Yscale','log')
grid on
axis([0.5 4.5 6e-4 1e-2])
ylabel('dRPE/dt, W/m^2')
set(gca,'XTickLabel',...
{'MPAS-O 120km','MPAS-O 60km','POP 1 deg','MOM 1 deg'})
group_legend={...
'\kappa_v=1e-5','\kappa_v=1e-6','\kappa_v=1e-7','\kappa_v=0',...
}
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.5 4 ])
fig=[char(dir) abc(1:nsims) '_DrpePertDtPOPBar'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
