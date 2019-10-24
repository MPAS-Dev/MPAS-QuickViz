% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:201];
dir={'c16'}; abc='abcd'; dims=1*[4 58];
grid_spacing=5e3;
nu_h=[.01 1 15 150];
% load('data/m52_igw_m52abcd.mat')

time_fields=[1:201];
dir={'m52'}; abc='abcdefghijklmnopqrstuvwx'; dims=1*[2 57];
grid_spacing=5e3;
nu_h=[.01 1 15 150];
load('data/m52_igw_m52a-x.mat') 

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'internal wave test',...
	  };

min_n=2;
max_n =11;

if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt = zeros(1,length(abc));
keMeanTime = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt(j),keMeanTime(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt ' char(dir) abc(j) ': %e \n'],meanDrpeDt(j));
end
end

% save(['data/' char(dir) abc '_rpe.mat'])

cl=[...
  0    0    0   ; ... % black
  0    0.6  0.9 ; ... % blue
  1    0    0   ; ... % red
  0    0.4  0.1 ; ... % dark green
  0.1  1    0.1 ; ... % green
  1.0  0.84 0   ; ... % gold  
  1    0.4  0   ; ... % orange
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];
blue = [0    0.6  0.9];

%%%%%%%%%%%%%%%%%%%%%
group_legend={'MPAS-O z-star','MITGCM','MOM',... 
	       'MPAS-O \tau_{Dlf}=0.1',...
	       'MPAS-O \tau_{Dlf}=1',...
	       'MPAS-O \tau_{Dlf}=10',...
	       'MPAS-O \tau_{Dlf}=100',...
	       'MPAS-O \tau_{Dlf}=1000',...
	      };
tau_Dlf = [.1 1 5 10 100 1000];

% for dam break, velocity scale is from eqn 5:
% or:
%rho1=999;
%rho2=997;
%vel_scale = 1/2*sqrt(9.8*500*(rho1-rho2)/((rho1+rho2)/2)) 



%%%%%% plot dRPE/dt versus grid Reynolds number
figure(10); clf
m=length(nu_h); % # experiments in set
marker='s*';
for j=2:2
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  vel_scale = sqrt(2*keMeanTime(iBeg:iEnd));
  gridRe = grid_spacing*vel_scale./nu_h;
  h=loglog(gridRe,meanDrpeDt(iBeg:iEnd),'-*k');
  set(h,'Color','k','LineWidth',1)
  set(h,'Marker',marker(j))
  hold on
end
grid on
axis([1 1e5 3e-5 1e-3])

load '~/mpas/email/130911_Ilicak_fig_data/fig14.mat'
h=loglog(mitgcm_Re,mitgcm_drpe_dt,'^g');
set(h,'MarkerFaceColor','g','Color','g');
h=loglog(mom_Re,mom_drpe_dt,'ob');
set(h,'MarkerFaceColor','b','Color','b');

if (1==1)
for j=3:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  vel_scale = sqrt(2*keMeanTime(iBeg:iEnd));
  gridRe = grid_spacing*vel_scale./nu_h;
  h=loglog(gridRe,meanDrpeDt(iBeg:iEnd),'-*');
  set(h,'Color',cl(j-1,:),'LineWidth',1)
  hold on
end
end

xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
legend(group_legend, 'location','NorthEastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 7 4])
fig=[char(dir) abc '_DrpeDt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, many viscosities

cl=[...
  0    0    0   ; ... % black
  1    0    0   ; ... % red
  0    0.6  0.9 ; ... % blue
  0.1  1    0.1 ; ... % green
  1.0  0.84 0   ; ... % gold  
  1    0.4  0   ; ... % orange
  0    0.4  0.1 ; ... % dark green
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];

group_legend={'\nu_h=0.01, Re_\Delta=5e4','\nu_h=1, Re_\Delta=500','\nu_h=15, Re_\Delta=30','\nu_h=150, Re_\Delta=2'};
figure(11); clf
nsims=4;
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot(time/3600/24,rpeNorm(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
end
grid on
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
axis tight
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_rpe_nu'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, one only

group_legend={'MPAS-O z-star','MITGCM','MOM','MPAS-O isopycnal.'...
	      };
figure(12); clf
nsims=1;
for j=1:nsims
  h=plot(time/3600/24,rpeNorm(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
end

load '~/mpas/email/130911_Ilicak_fig_data/fig13e.mat'
h=plot(mitgcm_time_days,mitgcm_rpe,'r-');
set(h,'Color','r','LineWidth',1);
h=plot(mom_time_days,mom_rpe,'-');
set(h,'Color',blue,'LineWidth',1);
h=plot(gold_time_days,gold_rpe,'--');
set(h,'Color',[ 0    0.4  0.1],'LineWidth',1); % dark yellow

grid on
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')
axis([0 200 -1e-6 9e-6])
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir) abc]);
fig=[char(dir) abc(1:nsims) '_rpe1'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);


return

