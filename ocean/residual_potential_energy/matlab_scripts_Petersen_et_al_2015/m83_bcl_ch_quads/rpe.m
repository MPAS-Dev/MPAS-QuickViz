% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

load '~/mpas/email/130911_Ilicak_fig_data/fig17.mat'

cl=[...
  0    0    0   ; ... % black
  1    0    0   ; ... % red
  0    0.6  0.9 ; ... % blue
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
marker='*s^v';
title_txt='bc';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 km

dir={'m81'}; abc='qmnop'; dims=1*[160 576];
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
else
  load('/data/pop_1km_140127.mat')
end

%meanDrpeDt1km100day = (rpe1km2-rpe1km)/320/86400;

% computed from a different file, see m62_bcl_ch_paper/rpe.m
% meanDrpeDt1km100day m62q: 1.593383e-04
% m62q keMeanTime1km =  0.001108653774138
meanDrpeDt1km100day(1) = 1.593383e-04;
keMeanTime1km(1) = 0.001108653774138;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 km 320 day

dir={'m83'}; abc='ghijkGHIJK'; dims=1*[40 143]; % change after m51k is done
nu_h=[1 5 10 20 200];
grid_spacing=4e3;
time_fields=[1:33];min_n=1;max_n = 33;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDtQuadHex4km = zeros(1,length(abc));
keMeanTimeQuadHex4km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDtQuadHex4km(j),keMeanTimeQuadHex4km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDtQuadHex4km ' char(dir) abc(j) ': %e \n'],meanDrpeDtQuadHex4km(j));
end
  save('data/m83QuadHex4km.mat','meanDrpeDtQuadHex4km','keMeanTimeQuadHex4km','dir','abc','nu_h')
else
  load('data/m83QuadHex4km.mat')
end

%%%%%% plot dRPE/dt versus grid Reynolds number

m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTimeQuadHex4km);
blue=[ 0    0.6  0.9];
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDtQuadHex4km(iBeg:iEnd),'k');
  set(h,'Color',   [ 0    0.6  0.9]) ;% blue
  set(h,'Marker',marker(j+2),'LineWidth',1)
  set(h,'LineStyle','none')
  hold on
end

h=loglog(mitgcm_4km_Re,mitgcm_4km_drpe_dt,'^b');
set(h,'MarkerFaceColor','g','Color','g');
h=loglog(mom_4km_Re,mom_4km_drpe_dt,'om');
set(h,'MarkerFaceColor','b','Color','b');

grid on
axis([.7 1e3 1e-6 3e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
group_legend={'MPAS-O z quad','MPAS-O hex 2 ord',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])


return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10 km hex

time_fields=[1:33];min_n=1;max_n = 33;

dir={'m83'}; abc='abcdeABCDE'; dims=1*[16 56];

nu_h=[1 5 10 20 200];
grid_spacing=10e3;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'baroclinic channel',...
	  };
if (1==1)
rpeNorm = zeros(length(time_fields),length(abc));
meanDrpeDtQuadHex10km = zeros(1,length(abc));
keMeanTimeQuadHex10km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
  [time,rpeTot,rpeNorm(:,j),DrpeDt,meanDrpeDtQuadHex10km(j),keMeanTimeQuadHex10km(j),vertTransportMean(j),vertTransportZ(:,j)] ...
    = sub_rpe(wd,dir,abc(j),netcdf_file, ...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDtQuadHex10km ' char(dir) abc(j) ': %e \n'],meanDrpeDtQuadHex10km(j));
end
  save('data/m83QuadHex10km.mat','meanDrpeDtQuadHex10km','keMeanTimeQuadHex10km','dir','abc','nu_h')
else
  load('data/m83QuadHex10km.mat')
end

%%%%%% plot dRPE/dt versus grid Reynolds number

figure(1)
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTimeQuadHex10km);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDtQuadHex10km(iBeg:iEnd),'k');
  set(h,'Marker',marker(j+2),'LineWidth',1)
  set(h,'LineStyle','none');
  hold on
end

h=loglog(mitgcm_10km_Re,mitgcm_10km_drpe_dt,'^b');
set(h,'MarkerFaceColor','g','Color','g');
h=loglog(mom_10km_Re,mom_10km_drpe_dt,'om');
set(h,'MarkerFaceColor','b','Color','b');

%title('MPAS-O Baroclinic channel, 10km')
grid on
axis([.7 1e3 1e-6 3e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
group_legend={'MPAS-O z* 100 day','MPAS-O z 100 day',...
	      'MPAS-O z* 320 day','MPAS-O z 320 day',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])


