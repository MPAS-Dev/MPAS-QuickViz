% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:40];
dir={'m66'}; abc='abcdef'; dims=1*[4 230];
min_n = 14;
max_n = 38;
group_legend={'MPAS-O, RPE no topo','MITGCM','MOM'};

time_fields=[1:40];
dir={'m66'}; abc='abcdefghijklABCDEFmnopqrstuvwx'; dims=1*[4 230];
min_n = 14;
max_n = 38;
group_legend={'MPAS-O, z-full cell','MPAS-O, z-pbc','MPAS-O, \sigma','MITGCM','MOM'};

time_fields=[1:40];
dir={'m66'}; abc='klEF'; dims=1*[4 230];
min_n = 14;
max_n = 38;
group_legend={'MPAS-O, z-full cell','MPAS-O, z-pbc','MPAS-O, \sigma','MITGCM','MOM'};

time_fields=[1:40];
dir={'c03'}; abc='lm'; dims=1*[4 230];
min_n = 1 ;
max_n = 7 ;
group_legend={'MPAS-O, z-full cell','MPAS-O, z-pbc','MPAS-O, \sigma','MITGCM','MOM'};


netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'overflow',...
	  };

%load('data/m66_overflow_131021.mat')
if (1==1)
load '~/mpas/email/130911_Ilicak_fig_data/fig10.mat'

rpeNorm = zeros(length(time_fields),length(abc));
headLocation = zeros(length(time_fields),length(abc));
meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 100;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt(j),keMeanTime(j),vertTransportMean(j),vertTransportZ(:,j),...
  	headLocation(:,j),DrpeDx ] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt ' char(dir) abc(j) ': %e \n'],meanDrpeDt(j));
end
save([char(dir) abc '_dRPEdx.mat'])
end

cl=[...
  0    0    0   ; ... % black
  1    0    0   ; ... % red
  0    0.6  0.9 ; ... % blue
  0.1  1    0.1 ; ... % green
  1    0.4  0   ; ... % orange
  0    0.4  0.1 ; ... % dark green
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
    
    
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];

%%%%%%%%%%%%%%%%%%%%%
%nu_h=[10.^[-2:2] 200];
nu_h=[10.^[3]];
%nu_h=[10.^[-2:3]];
tau_Dlf = [.1 1 5 10 100 1000];
grid_spacing=1000;

%vel_scale = sqrt(2*keMeanTime);
% for overflow, velocity scale is from eqn 5:
% or:
rho1=999;
rho2=997;
vel_scale = 1/2*sqrt(9.8*500*(rho1-rho2)/((rho1+rho2)/2)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE versus X

group_legend={'z* \nu_h=100','z* \nu_h=1000','\sigma \nu_h=100','\sigma \nu_h=1000',...
	      };
figure(10); clf
for j=1:2
  h=plot(headLocation(:,j)/1e3,rpeNorm(:,j),'-o');
  set(h,'Color',cl(2,:),'LineWidth',1)
  if j==1; set(h,'MarkerFaceColor',cl(2,:)); end
  hold on
end

hold on
%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:15];
dir={'c03'}; abc='rs'; dims=1*[4 230];
min_n = 1 ;
max_n = 7 ;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'overflow',...
	  };


%load('data/m66_overflow_131021.mat')
if (1==1)
load '~/mpas/email/130911_Ilicak_fig_data/fig10.mat'

rpeNorm = zeros(length(time_fields),length(abc));
headLocation = zeros(length(time_fields),length(abc));
meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 100;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt(j),keMeanTime(j),vertTransportMean(j),vertTransportZ(:,j),...
  	headLocation(:,j),DrpeDx ] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt ' char(dir) abc(j) ': %e \n'],meanDrpeDt(j));
end
save([char(dir) abc '_dRPEdx.mat'])
end


%%%%%%%%%%%%%%%%%%%%%
%nu_h=[10.^[-2:2] 200];
nu_h=[10.^[3]];
%nu_h=[10.^[-2:3]];
tau_Dlf = [.1 1 5 10 100 1000];
grid_spacing=1000;

%vel_scale = sqrt(2*keMeanTime);
% for overflow, velocity scale is from eqn 5:
% or:
rho1=999;
rho2=997;
vel_scale = 1/2*sqrt(9.8*500*(rho1-rho2)/((rho1+rho2)/2)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE versus X

for j=1:2
  h=plot(headLocation(:,j)/1e3,rpeNorm(:,j),'-o');
  set(h,'Color',cl(3,:),'LineWidth',1)
  if j==1; set(h,'MarkerFaceColor',cl(3,:)); end
end

grid on
xlabel('head location, km')
ylabel('(RPE-RPE(0))/RPE(0)')
axis([0 180 -.5e-6 7e-6])
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

fig=[char(dir) abc '_DrpeDx'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);




return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% head location vs time

group_legend={'z-pbc, \nu_h=100','z-pbc, \nu_h=1000','\sigma, \nu_h=100','\sigma, \nu_h=1000',...
	      };
figure(11); clf
for j=1:2
  h=plot(time/3600,headLocation(:,j)/1e3,'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
end

for j=3:4
  h=plot(time/3600,headLocation(:,j)/1e3,'--');
  set(h,'Color',cl(j-2,:),'LineWidth',1)
  hold on
end

grid on
xlabel('time, hours')
ylabel('head location, km')
legend(group_legend,'location','SouthEast')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

fig=[char(dir) abc '_head_location'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
return


%%%%%% plot dRPE/dt versus grid Reynolds number

figure(2); clf
m=length(nu_h); % # experiments in set
marker='s*';
for j=1:1
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale./nu_h
  h=loglog(gridRe,meanDrpeDt(iBeg:iEnd),'-*k');
  set(h,'Color',cl(j,:))
  hold on
end
grid on
axis([1 1e6 1e-1 3e1])

h=loglog(mitgcm_Re,mitgcm_drpe_dt,'^k');
set(h,'MarkerFaceColor','k');
%set(h,'Color',   [ 1    0    0],'LineWidth',1) ;% red
h=loglog(mom_Re,mom_drpe_dt,'ok');
set(h,'MarkerFaceColor','k');

for j=4:5
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale./nu_h
  h=loglog(gridRe,meanDrpeDt(iBeg:iEnd),'s');
  set(h,'Color',cl(j-3,:))
  hold on
end


xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
legend(group_legend,...
       'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4 3])
fig=[char(dir) abc '_DrpeDt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
%fig=[char(dir) abc '_DrpeDt'];

% for printed page:
%set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
%	'PaperPosition',[1 3 6.5 5.0])
%subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir)]);% abc]);
%unix(['convert f/' fig '.eps f/' fig '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, many viscosities

group_legend={'\nu_h=0.01','\nu_h=0.1','\nu_h=1','\nu_h=10','\nu_h=100','\nu_h=1000',...
	      };
figure(10); clf
for j=1:6
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot(time/3600,rpeNorm(:,j+6),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
end
grid on
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
axis tight
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

fig=[char(dir) abc(7:12) '_rpe_nu'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, many viscosities

group_legend={'z*-full cell','z*-pbc','sigma','z-full cell','z-pbc'...
	      };
figure(20); clf
cases=4:6:30;
abc(cases)
cli=1
for j=cases
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot(time/3600,rpeNorm(:,j),'-');
  set(h,'Color',cl(cli,:),'LineWidth',1)
  cli=cli+1;
  hold on
end
grid on
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
axis tight
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

fig=[char(dir) abc(cases) '_rpe_nu'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);



return
%%%%%%%%%%%%% add non-Rich to plot

time_fields=[1:40];
dir={'m57'}; abc='d'; dims=1*[2 230];
min_n =  7;
max_n = 11;

if (1==1)

rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 100;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt(j),keMeanTime(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt ' char(dir) abc(j) ': %e \n'],meanDrpeDt(j));
end
end

figure(11);
nsims=1;
for j=1:nsims
  h=plot(time/3600,rpeNorm(:,j),'--');
  set(h,'Color',cl(4,:),'LineWidth',1)
  hold on
end

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

fig=[char(dir) abc(1:nsims) '_rpe_nu_w_nonRich'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
