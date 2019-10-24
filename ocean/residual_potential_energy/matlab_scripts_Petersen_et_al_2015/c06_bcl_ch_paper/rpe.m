% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

load '~/mpas/email/130911_Ilicak_fig_data/fig17.mat'

title_txt='bc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 km 320 day

%dir={'m62'}; abc='ghijkwxyzZ'; dims=[40 143];
%dir={'c07'}; abc='abcdefghijklmnoFGHIJ'; dims=[40 143];
%dir={'c15'}; abc='klmn'; dims=[40 143]; 
dir={'c12'}; abc='ABCDEFGHIJpqrstuvwxy'; dims=[40 143]; 
dir={'c15'}; abc='FGHIJPQRSTabcdefghijABCDE'; dims=[40 143]; 
nu_h=[1 5 10 20 200];
grid_spacing=4e3;
time_fields=[1:32];min_n=1;max_n = 32;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt4km320day = zeros(1,length(abc));
keMeanTime4km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt4km320day(j),keMeanTime4km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt4km320day ' char(dir) abc(j) ': %e \n'],meanDrpeDt4km320day(j));
end
else
  load(['data/' char(dir) abc '_rpe.mat'])
end
% save([char(dir) abc '_rpe.mat'])

cl=[...
  0    0.4  0.1 ; ... % dark green
  1    0    0   ; ... % red
  0    0.6  0.9 ; ... % blue
  0.1  1    0.1 ; ... % green
  0    0    0   ; ... % black
  1    0.4  0   ; ... % orange
  0.6  0.6  0.6 ; ... % grey
  1    0.8  0   ; ... % yellow
  0.6  0    0.8 ; ... % purple
  1    0    1   ; ... % pink
  0    1    1   ; ... % aquablue
  0.4  0.4  0   ; ... % army green
];
marker='sxo+*';

%%%%%% plot dRPE/dt versus grid Reynolds number

m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime4km);
blue=[ 0    0.6  0.9];
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt4km320day(iBeg:iEnd),'k');
  set(h,'Color', cl(j,:)) ;% blue
  set(h,'Marker',marker(j),'LineWidth',1)
  %if j>1; set(h,'LineStyle','--'); end
  hold on
end


%%%%%%%%%%%%% plot POP data

load('data/pop_4km_140127.mat')

vel_scale = sqrt(2*keMeanTime4kmPOP);
nu_h=[1 5 10 20 200];
  iBeg=1;
  iEnd=5;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt4kmPOP(iBeg:iEnd),'rs');
set(h,'MarkerFaceColor','r','Color','r');

  
h=loglog(mitgcm_4km_Re,mitgcm_4km_drpe_dt,'^k');
set(h,'MarkerFaceColor','g','Color','g');
h=loglog(mom_4km_Re,mom_4km_drpe_dt,'ok');
set(h,'MarkerFaceColor','b','Color','b');

grid on
%axis([.1 2e3 1e-6 3e-3])
axis([.7 1e3 2e-7 2e-2])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
h=text(1,1.3e-2,'4 km resolution');
set(h,'FontSize',12,'FontWeight','bold')
group_legend={'MPAS 1st order advection','MPAS FCT: F_L=1st F_H=25% 1st, 75% 2nd',...
	      'MPAS FCT: F_L=1st F_H=5% 1st, 95% 2nd','MPAS FCT: F_L=1st F_H=2nd',...
	      'MPAS standard FCT: F_L=1st F_H=3rd',...
	      'POP               MITGCM               MOM',
	     };
legend(group_legend,'location','SouthEast')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5 4.5])
%	'PaperPosition',[0 0 5 3.3])

fig=[char(dir) '_DrpeDt_4km_order2blend_part_legend'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10 km 320 days

time_fields=[1:33];min_n=1;max_n = 33;

%dir={'m62'}; abc='abcderstuv'; dims=[16 56];
dir={'c06'}; abc='abcdeghijk'; dims=[16 58];
dir={'c06'}; abc='mnopq'; dims=[16 56];

nu_h=[1 5 10 20 200];
grid_spacing=10e3;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'baroclinic channel',...
	  };
if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt10km320day = zeros(1,length(abc));
keMeanTime10km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt10km320day(j),keMeanTime10km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt10km320day ' char(dir) abc(j) ': %e \n'],meanDrpeDt10km320day(j));
end
end

%%%%%% plot dRPE/dt versus grid Reynolds number

figure(1)
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime10km);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt10km320day(iBeg:iEnd),'k');
  set(h,'Marker',marker(j),'LineWidth',1)
  %set(h,'LineStyle','none');
  hold on
end

h=loglog(mitgcm_10km_Re,mitgcm_10km_drpe_dt,'^g');
set(h,'MarkerFaceColor','g');
%set(h,'Color',   [ 1    0    0],'LineWidth',1) ;% red
h=loglog(mom_10km_Re,mom_10km_drpe_dt,'ob');
set(h,'MarkerFaceColor','b');

%title('MPAS-O Baroclinic channel, 10km')
grid on
axis([.1 2e3 1e-6 3e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
group_legend={'MPAS-O z* 100 day','MPAS-O z 100 day',...
	      'MPAS-O z* 320 day','MPAS-O z 320 day',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 km, m62q only

dir={'m62'}; abc='q'; dims=[160 576];
nu_h=[1];
grid_spacing=1e3;

netcdf_file = 'output.0000-01-01_00.00.00_day1-320.nc';
time_fields=[1:2];min_n=1;max_n=2;

if (1==1)
  clear time DrpeDt rpeTot
rpe1km = zeros(1,length(abc));
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt1km100day = zeros(1,length(abc));
keMeanTime1km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time, rpeTot, rpeNorm(:,j), DrpeDt,meanDrpeDt1km100day(j),  keMeanTime1km(j),vertTransportMean(j),vertTransportZ(:,j)]...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt1km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt1km100day(j));
end
end

% meanDrpeDt1km100day m62q: 1.593383e-04
% m62q keMeanTime1km =  0.001108653774138

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 km

dir={'m62'}; abc='qmnop'; dims=[160 576];
nu_h=[1 5 10 20 200];
grid_spacing=1e3;

netcdf_file = 'output.0000-01-01_00.00.00.nc';
time_fields=[1:11];min_n=1;max_n=11;

if (1==2)
  clear time DrpeDt rpeTot
rpe1km = zeros(1,length(abc));
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt1km100day = zeros(1,length(abc));
keMeanTime1km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time, rpeTot, rpeNorm(:,j), DrpeDt,meanDrpeDt1km100day(j),  keMeanTime1km(j),vertTransportMean(j),vertTransportZ(:,j)]...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt1km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt1km100day(j));
end
end

%meanDrpeDt1km100day = (rpe1km2-rpe1km)/320/86400;

%%%%%% plot dRPE/dt versus grid Reynolds number

h=loglog(mitgcm_1km_Re,mitgcm_1km_drpe_dt,'^k');
set(h,'MarkerFaceColor','r','Color','r');
h=loglog(mom_1km_Re,mom_1km_drpe_dt,'ok');
set(h,'MarkerFaceColor','r','Color','r');

m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime1km);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt1km100day(iBeg:iEnd),'r');
  %set(h,'Color',cl(j,:))
  set(h,'Marker',marker(j),'LineWidth',1)
  if j==1; set(h,'LineStyle','-'); end
  hold on
end

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5 3.3])

fig=[char(dir) '_DrpeDt_small'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 6 4])

fig=[char(dir) '_DrpeDt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10 km 100 day

time_fields=[1:11];min_n=1;max_n = 11;

dir={'m62'}; abc='abcderstuv'; dims=[16 56];

nu_h=[1 5 10 20 200];
grid_spacing=10e3;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'baroclinic channel',...
	  };
if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt10km100day = zeros(1,length(abc));
keMeanTime10km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt10km100day(j),keMeanTime10km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt10km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt10km100day(j));
end
end

%%%%%% plot dRPE/dt versus grid Reynolds number

figure(1); clf;
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime10km);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt10km100day(iBeg:iEnd),'k');
  set(h,'Marker',marker(j),'LineWidth',1)
  if j>1; set(h,'LineStyle','none'); end
  hold on
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10 km 320 days

time_fields=[1:33];min_n=1;max_n = 33;

dir={'m62'}; abc='abcderstuv'; dims=[16 56];

nu_h=[1 5 10 20 200];
grid_spacing=10e3;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'baroclinic channel',...
	  };
if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt10km320day = zeros(1,length(abc));
keMeanTime10km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt10km320day(j),keMeanTime10km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt10km320day ' char(dir) abc(j) ': %e \n'],meanDrpeDt10km320day(j));
end
end

%%%%%% plot dRPE/dt versus grid Reynolds number

figure(1)
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime10km);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt10km320day(iBeg:iEnd),'k');
  set(h,'Marker',marker(j+2),'LineWidth',1)
  set(h,'LineStyle','none');
  hold on
end

h=loglog(mitgcm_10km_Re,mitgcm_10km_drpe_dt,'^k');
set(h,'MarkerFaceColor','k');
%set(h,'Color',   [ 1    0    0],'LineWidth',1) ;% red
h=loglog(mom_10km_Re,mom_10km_drpe_dt,'ok');
set(h,'MarkerFaceColor','k');

%title('MPAS-O Baroclinic channel, 10km')
grid on
axis([.1 2e3 1e-6 3e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
group_legend={'MPAS-O z* 100 day','MPAS-O z 100 day',...
	      'MPAS-O z* 320 day','MPAS-O z 320 day',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 km

dir={'m62'}; abc='ghijkwxyzZ'; dims=[40 143]; % change after m51k is done
nu_h=[1 5 10 20 200];
grid_spacing=4e3;
time_fields=[1:11];min_n=1;max_n = 11;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDt4km100day = zeros(1,length(abc));
keMeanTime4km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDt4km100day(j),keMeanTime4km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt4km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt4km100day(j));
end
end

%%%%%% plot dRPE/dt versus grid Reynolds number

m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime4km);
blue=[ 0    0.6  0.9];
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt4km100day(iBeg:iEnd),'k');
  set(h,'Color',   [ 0    0.6  0.9]) ;% blue
  set(h,'Marker',marker(j),'LineWidth',1)
  if j>1; set(h,'LineStyle','none'); end
  hold on
end

