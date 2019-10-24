% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';
load '~/mpas/email/130911_Ilicak_fig_data/fig20.mat'

read_30km_data= 0;
read_15km_data= 0;
min_n=1; max_n=1;

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
%       'meanDensity','varDensity','rpeMean','rpePert');
end

group_legend={...
'MPAS-O 30km \kappa_v=1e-5','MPAS-O 30km \kappa_v=1e-6','MPAS-O 30km \kappa_v=1e-7','MPAS-O 30km \kappa_v=0',...
'MOM 1/4^o \kappa_v=1e-5','MOM 1/4^o \kappa_v=1e-6','MOM 1/4^o \kappa_v=0',...
	      };
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=0.1',...
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=1',...
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=5',...
%'MPAS-O 30km \kappa_v=0 \tau_{Dlf}=10',...

cl=[...
  0    0.6  0.9 ; ... % blue
  1    0    0   ; ... % red
  0.1  1    0.1 ; ... % green
  0    0    0   ; ... % black
  1    0.4  0   ; ... % orange
  0    0.4  0.1 ; ... % dark green
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];

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
       'meanDensity','varDensity','rpeMean','rpePert');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 120 km

time_fields=[1:71];
%dir={'m55'}; abc='bcdeqrs'; nVertLevels = 40;
dir={'m55'}; abc='a'; nVertLevels = 40;

netcdf_file = 'output.0010-01-01_00.00.00_rpe_vars.nc';

group_legend={...
'MPAS-O 120km \kappa_v=0','MPAS-O 120km \kappa_v=1e-7','MPAS-O 120km \kappa_v=1e-6','MPAS-O 120km \kappa_v=1e-5',...
'MPAS-O 120km \kappa_v=0 \tau_{Dlf}=1',...
'MPAS-O 120km \kappa_v=0 \tau_{Dlf}=5',...
'MPAS-O 120km \kappa_v=0 \tau_{Dlf}=10',...
	      };

group_legend={...
'MPAS-O 120km forced'
	      };

title_txt={'global spin-down' };

if (read_120km_data)

  rpeNorm = zeros(length(time_fields),length(abc));
  DrpeDt = zeros(length(time_fields),length(abc));
  keMeanVolume = zeros(length(time_fields),length(abc));
  vertTransportVolume = zeros(length(time_fields),length(abc));
  vertTransportZ = zeros(nVertLevels,length(abc));
  for j=1:length(abc)
    [time,rpeTot,rpeNorm(:,j),DrpeDt(:,j),keMeanVolume(:,j),vertTransportVolume(:,j),vertTransportMeanZ(:,j)] ...
      = sub_rpe(wd,dir,abc(j),netcdf_file, ...
      time_fields,min_n,max_n,char(title_txt));
  end
  save(['data/' char(dir) abc '_global_spindown.mat'])
else
  load(['data/' char(dir) abc '_global_spindown.mat'])
end

meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
for j=1:length(abc)
  meanDrpeDt(j) = mean(DrpeDt(min_n:max_n));
  vertTransportMean(j) = mean(vertTransportVolume(min_n:max_n));
end



