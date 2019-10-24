% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';
%load '~/mpas/email/130911_Ilicak_fig_data/fig20.mat'

read_120km_data=0;
read_60km_data=0;
read_30km_data= 0;
read_15km_data= 0;
min_n=1; max_n=1;
nrez=4;

cl=[...
  0    0    0   ; ... % black
  1    0    0   ; ... % red
  0    0.6  0.9 ; ... % blue
  0.05  .8    0.1 ; ... % mid-green
  0.1  1    0.1 ; ... % green
  0    0.4  0.1 ; ... % dark green
  1    0.4  0   ; ... % orange
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];

meanDrpePertDtAll = zeros(nrez,3);
gridReAll = zeros(nrez,3);
nu_h_del4All = zeros(nrez,3);

% compute grid Re using hyperviscosity:
% http://physics.stackexchange.com/questions/21108/reynolds-number-with-hyper-viscosity
% Re = U L^{2H-1}/nu_h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 120 km

time_fields=[2:8];
dir={'m89'}; abc='bcaed'; nVertLevels = 40; rezIndex=1;

netcdf_file = 'output.0001-01-11_00.00.00.nc';
nu_h_del4=2.6e13*[8 sqrt(8) 1 1/sqrt(8) 1/8];
grid_spacing = 120e3;

title_txt={'global spin-down' };

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  multiple_files=false;

if (read_120km_data)

  for j=1:length(abc)
    [time,rpeTot(:,j),rpeNorm(:,j),DrpeDt(:,j),keMeanVolume(:,j), ... %vertTransportVolume(:,j),vertTransportMeanZ(:,j), ...
	  rpeMean(:,j),rpePert(:,j)] ...
	= sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt),multiple_files);
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'],'time','DrpeDt','rpeTot','rpeNorm',...
       'rpeMean','rpePert', 'keMeanVolume');
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
  meanDrpePertDtAll(rezIndex,j) = mean(DrpePertDt(:,j));

  keMeanTime(j) = mean(keMeanVolume(:,j));
  vel_scale(j) = sqrt(2*keMeanTime(j));
  gridReAll(rezIndex,j) = grid_spacing^3*vel_scale(j)/nu_h_del4(j);
  nu_h_del4All(rezIndex,j) = nu_h_del4(j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 60 km
time_fields=[2:8];
dir={'m89'}; abc='ghfji'; nVertLevels = 40;rezIndex=2;
netcdf_file = 'output.0001-01-11_00.00.00.nc';
nu_h_del4=3.2e12*[8 sqrt(8) 1 1/sqrt(8) 1/8];
grid_spacing = 60e3;

title_txt={'global spin-down' };

if (read_60km_data)

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  multiple_files=false;
  for j=1:length(abc)
    [time,rpeTot(:,j),rpeNorm(:,j),DrpeDt(:,j),keMeanVolume(:,j), ... %vertTransportVolume(:,j),vertTransportMeanZ(:,j), ...
	  rpeMean(:,j),rpePert(:,j)] ...
	= sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt),multiple_files);
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'],'time','DrpeDt','rpeTot','rpeNorm',...
       'rpeMean','rpePert', 'keMeanVolume');
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
  meanDrpePertDtAll(rezIndex,j) = mean(DrpePertDt(:,j));

  keMeanTime(j) = mean(keMeanVolume(:,j));
  vel_scale(j) = sqrt(2*keMeanTime(j));
  gridReAll(rezIndex,j) = grid_spacing^3*vel_scale(j)/nu_h_del4(j);
  nu_h_del4All(rezIndex,j) = nu_h_del4(j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 30 km
time_fields=[1:6];
dir={'m89'}; abc='lmkon'; nVertLevels = 40;rezIndex=3;

netcdf_file = 'output.0001-01-11_00.00.00_rpe_vars';
nu_h_del4=4.0e11*[8 sqrt(8) 1 1/sqrt(8) 1/8];
grid_spacing = 30e3;

if (read_30km_data)

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  multiple_files=true;
  for j=1:length(abc)
    [time,rpeTot(:,j),rpeNorm(:,j),DrpeDt(:,j),keMeanVolume(:,j), ... %vertTransportVolume(:,j),vertTransportMeanZ(:,j), ...
	  rpeMean(:,j),rpePert(:,j)] ...
	= sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt),multiple_files);
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'],'time','DrpeDt','rpeTot','rpeNorm',...
       'rpeMean','rpePert', 'keMeanVolume');
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
  meanDrpePertDtAll(rezIndex,j) = mean(DrpePertDt(:,j));

  keMeanTime(j) = mean(keMeanVolume(:,j));
  vel_scale(j) = sqrt(2*keMeanTime(j));
  gridReAll(rezIndex,j) = grid_spacing^3*vel_scale(j)/nu_h_del4(j);
  nu_h_del4All(rezIndex,j) = nu_h_del4(j);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 15 km
time_fields=[1:6];
dir={'m89'}; abc='rps'; nVertLevels = 40;rezIndex=4;

netcdf_file = 'output.0001-01-11_00.00.00_rpe_vars';
nu_h_del4=5.0e10*[sqrt(8) 1 1/8];
grid_spacing = 15e3;

if (read_15km_data)

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  multiple_files=true;
  for j=1:length(abc)
    [time,rpeTot(:,j),rpeNorm(:,j),DrpeDt(:,j),keMeanVolume(:,j), ... %vertTransportVolume(:,j),vertTransportMeanZ(:,j), ...
	  rpeMean(:,j),rpePert(:,j)] ...
	= sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt),multiple_files);
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'],'time','DrpeDt','rpeTot','rpeNorm',...
       'rpeMean','rpePert', 'keMeanVolume');
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
  meanDrpePertDtAll(rezIndex,j) = mean(DrpePertDt(:,j));

  keMeanTime(j) = mean(keMeanVolume(:,j));
  vel_scale(j) = sqrt(2*keMeanTime(j));
  gridReAll(rezIndex,j) = grid_spacing^3*vel_scale(j)/nu_h_del4(j);
  nu_h_del4All(rezIndex,j) = nu_h_del4(j);
end

%%%%%%%%%%%% dRPE/dt vs grid Re
% 141124: divide grid Reynolds number by a factor of 8
gridReAll = gridReAll/8.0;

figure(1)
h = loglog(gridReAll',meanDrpePertDtAll','-*');
hold on
grid on
legend('120km','60km','30km','15km','Location','NorthWest')

for j=1:4
  set(h(j),'Color',cl(j,:),'LineWidth',1)
end

% Add boxes to typical configuration:
kk=[3 3 3 2];
for j=1:4
  hb = loglog(gridReAll(j,kk(j)),meanDrpePertDtAll(j,kk(j)),'s');
  set(hb,'Color',cl(j,:),'MarkerSize',10)
end

xlabel('grid Reynolds Number')
ylabel('dRPE/dt, W/m^2')
axis([.05 15 4e-4 3e-3])
set(gca,'yTick',[1e-4*[1:9] 1e-3*[1 2 3]])

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
        'PaperPosition',[0 0 5 3.3])

fig=[char(dir) '_DrpeDt_vs_Re_tight_factor8'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return


%%%%%%%%%%%% bar chart

figure(2); clf;
h=bar(meanDrpePertDtAll,'ShowBaseLine', 'off')
set(h,'BaseValue',1e-4);
set(gca,'Yscale','log')
grid on
axis([0.5 4.5 1e-4 1e-2])
ylabel('dRPE/dt, W/m^2')
set(gca,'XTickLabel',...
{'MPAS-O 120km','MPAS-O 60km'})
group_legend={...
'\nu=high','\nu=med','\nu=low',...
}
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.5 4 ])
fig=[dir abc '_DrpePertDtBar'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);



return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 30 km

dir={'m53'}; abc='xyz'; nVertLevels = 40;
dir={'m53'}; abc='ps'; nVertLevels = 40;
time_fields=[0:9];
dir={'m53'}; abc='wvut'; nVertLevels = 40;
time_fields=[0:5:70]; 

netcdf_file = 'output.0010-01-01_00.00.00_rpe_vars';

title_txt={'global spin-down' };

if (read_30km_data)

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  for j=1:length(abc)
    [time,rpeTot(:,j),rpeNorm(:,j),DrpeDt(:,j),... %keMeanVolume(:,j),vertTransportVolume(:,j),vertTransportMeanZ(:,j), ...
	  rpeMean(:,j),rpePert(:,j)] ...
	= sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt));
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'])
%  ,'time','DrpeDt','rpeTot','rpeNorm',...
%       'rpeMean','rpePert');
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
end

meanDrpePertDtAll = zeros(4,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dRPE/dt with time

figure(28); clf
nsims=length(abc);
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=semilogy((time-time(1))/3600/24,DrpePertDt(:,j),'--');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
  meanDrpePertDtAll(1,j) = mean(DrpePertDt(5:length(time),j))
end
grid on
xlabel('time, days')
ylabel('dRPE/dt, W/m^2')
axis([0 70 6e-4 1e-2])
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 4])

if (1==2)
j=1;
h=semilogy(mom_quarter_kd_0_time,mom_quarter_kd_0_drpe_dt,'--');
set(h,'Color',cl(j,:),'LineWidth',1); j=j+2
h=semilogy(mom_quarter_kd_1e_6_time,mom_quarter_kd_1e_6_drpe_dt,'--');
set(h,'Color',cl(j,:),'LineWidth',1); j=j+1;
h=semilogy(mom_quarter_kd_1e_5_time,mom_quarter_kd_1e_5_drpe_dt,'--');
set(h,'Color',cl(j,:),'LineWidth',1); 

legend(group_legend,'location','EastOutside')

%subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_DrpePertdt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPEnorm with time

figure(23); clf
nsims=length(abc);
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot((time-time(1))/3600/24,rpeNorm(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',2)
  hold on
end
grid on
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
axis tight
legend(group_legend,'location','EastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 4])

subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_DrpeNormdt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time

figure(24); clf
nsims=length(abc);
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot((time-time(1))/3600/24,rpeTot(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',2)
  hold on
end
grid on
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
axis tight
legend(group_legend,'location','EastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 4])

subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_rpe_kappa'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 15 km

time_fields=[0:70];
dir={'m53'}; abc='srqp'; nVertLevels = 40;

netcdf_file = 'output.0010-01-01_00.00.00_rpe_vars';

title_txt={'global spin-down' };

if (read_15km_data)

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  rpeMean = zeros(length(time_fields),length(abc));
  rpePert = zeros(length(time_fields),length(abc));
  rpeTot = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  for j=4:length(abc)
    [time,rpeTot(:,j),rpeNorm(:,j),DrpeDt(:,j),... %keMeanVolume(:,j),vertTransportVolume(:,j),vertTransportMeanZ(:,j), ...
	  rpeMean(:,j),rpePert(:,j)] ...
	= sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt));
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'],'time','DrpeDt','rpeTot','rpeNorm',...
       'rpeMean','rpePert');
end

group_legend={...
'MPAS-O 30km \kappa_v=1e-5','MPAS-O 30km \kappa_v=1e-6','MPAS-O 30km \kappa_v=1e-7','MPAS-O 30km \kappa_v=0',...
'MOM 1/4^o \kappa_v=1e-5','MOM 1/4^o \kappa_v=1e-6','MOM 1/4^o \kappa_v=0',...
	      };
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=0.1',...
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=1',...
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=5',...
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=10',...

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
DrpePertDt = zeros(size(rpePert));
DrpeMeanDt = zeros(size(rpeMean));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
  DrpePertDt(:,j) = derivative_ord2(time,rpePert(:,j)');
  DrpeMeanDt(:,j) = derivative_ord2(time,rpeMean(:,j)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dRPE/dt with time

figure(28);
nsims=length(abc);
for j=1:3:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=semilogy((time-time(1))/3600/24,DrpePertDt(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
  meanDrpePertDtAll(2,j) = mean(DrpePertDt(20:length(time),j));
end

% fill in from
% /home/mpeterse/mpas/matlab/rpe/m53qr_global_spindown_15km
% ran by hand
meanDrpePertDtAll(2,2) = 0.9295e-3 % m53r, kappa=1e-6
meanDrpePertDtAll(2,3) = 0.7961e-3 % m53q, kappa=1e-7

grid on
xlabel('time, days')
ylabel('dRPE/dt, W/m^2')
axis([0 70 6e-4 1e-2])
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 4])

j=1;
h=semilogy(mom_quarter_kd_1e_5_time,mom_quarter_kd_1e_5_drpe_dt,'-.');
set(h,'Color',cl(j,:),'LineWidth',1); j=j+1
  meanDrpePertDtAll(4,1) = mean(mom_quarter_kd_1e_5_drpe_dt(20:69));
h=semilogy(mom_quarter_kd_1e_6_time,mom_quarter_kd_1e_6_drpe_dt,'-.');
set(h,'Color',cl(j,:),'LineWidth',1); j=j+2;
  meanDrpePertDtAll(4,2) = mean(mom_quarter_kd_1e_6_drpe_dt(20:69));
h=semilogy(mom_quarter_kd_0_time,mom_quarter_kd_0_drpe_dt,'-.');
set(h,'Color',cl(j,:),'LineWidth',1); 
  meanDrpePertDtAll(4,4) = mean(mom_quarter_kd_0_drpe_dt(20:69));

group_legend={...
'MPAS-O 30km \kappa_v=1e-5','MPAS-O 30km \kappa_v=1e-6','MPAS-O 30km \kappa_v=1e-7','MPAS-O 30km \kappa_v=0',...
'MPAS-O 15km \kappa_v=1e-5','MPAS-O 15km \kappa_v=0',...
'MOM 1/4^o \kappa_v=1e-5','MOM 1/4^o \kappa_v=1e-6','MOM 1/4^o \kappa_v=0',...
};
%'MPAS-O 15km \kappa_v=1e-5','MPAS-O 15km \kappa_v=1e-6','MPAS-O 15km \kappa_v=1e-7','MPAS-O 15km \kappa_v=0',...

legend(group_legend,'location','EastOutside')

%subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_drpePertdt'];
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
{'MPAS-O 30km','MPAS-O 15km','','MOM 1/4'})
group_legend={...
'\kappa_v=1e-5','\kappa_v=1e-6','\kappa_v=1e-7','\kappa_v=0',...
}
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.5 4 ])
fig=['m53srqp_DrpePertDtBar'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time

figure(23); clf
nsims=length(abc);
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot((time-time(1))/3600/24,rpeNorm(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',2)
  hold on
end
grid on
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
axis tight
legend(group_legend,'location','EastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 4])

subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_rpe_kappa'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dRPE/dt with time

figure(24); clf
nsims=length(abc);
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=semilogy((time-time(1))/3600/24,DrpeDt(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',2)
  hold on
end
grid on
xlabel('time, days')
ylabel('dRPE/dt, W/m^2')
axis([0 70 1e-4 3e-2])
legend(group_legend,'location','EastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 8 4])

%subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_drpedt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);




