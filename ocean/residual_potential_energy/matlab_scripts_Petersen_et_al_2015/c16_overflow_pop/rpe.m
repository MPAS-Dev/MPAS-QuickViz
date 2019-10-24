% script rpe
% compute resting potential energy, ala Ilicak ea 2012 for POP data
% Mark Petersen, LANL, Dec 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/'; % on ccscs5

load '~/mpas/email/130911_Ilicak_fig_data/fig17.mat'

%load('data/m62_bcl_ch_131028.mat')
%save(['data/' char(dir) abc '_dRPEdt.mat'])

cl=[...
  0    0    0   ; ... % black
  1    0    0   ; ... % red
  0    0.6  0.9 ; ... % blue
  0.1  1    0.1 ; ... % green
  1    0.4  0   ; ... % orange
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
marker='s*';
title_txt='bc';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 km POP

dir={'c16'}; abc='e'; 

nu_h=[1000];
grid_spacing=1e3;

hist_files = {...
't.overflow8a.24hours.nc'};

time_fields=[1:25]; % time in days
min_n=14;
max_n=25;

kmt_file = 'kmt.dzbc';
grid_file = 'grid';
depth_file = 'in_depths.dat';

title_txt={
    'overflow POP',...
	  };
if (1==2)
rpeNorm = zeros(length(time_fields),length(abc));
rpeTot = zeros(length(time_fields),length(abc));
keMeanVolume = zeros(length(time_fields),length(abc));
meanDrpeDt1km100day = zeros(1,length(abc));
keMeanTime1km = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 100;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
[time,rpeTot(:,j),rpeNorm(:,j),DrpeDt,meanDrpeDt1km(j),keMeanTime1km(j),keMeanVolume(:,j)] ...
    = sub_rpe(wd,dir,abc(j),hist_files, kmt_file, grid_file, depth_file, ...
    time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt1km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt1km100day(j));
end

end

%%%%%% plot RPE norm versus time

figure(2); clf;
m=length(nu_h); % # experiments in set
for j=1:length(abc)
  h=plot(time_fields,rpeNorm(:,j));
  set(h,'LineWidth',1,'Color',cl(j,:))
  if j>m; set(h,'LineStyle','--'); end
  grid on
  hold on
end
axis tight
title('POP z-level overflow')
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')

group_legend={'\nu_h=1000'};
legend(group_legend,'location','NorthWest')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])
fig=[char(dir) '_10km_RPEnorm'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);


return

%%%%%% plot dRPE/dt versus grid Reynolds number
%load('data/pop_1km_140110.mat')   

figure(11);
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime1km);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe(1:4),meanDrpeDt1km(1:4),'r');
  set(h,'Marker',marker(j),'LineWidth',1)
  if j~=2; set(h,'LineStyle','none'); end
  grid on
  hold on
end

h=loglog(mitgcm_1km_Re,mitgcm_1km_drpe_dt,'^b');
set(h,'MarkerFaceColor','g','Color','g');
h=loglog(mom_1km_Re,mom_1km_drpe_dt,'om');
set(h,'MarkerFaceColor','b','Color','b');

%title('MPAS-O Baroclinic channel, 1km')
grid on
axis([.7 1e3 1e-6 3e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
h=text(1,1.5e-3,'1 km resolution')
set(h,'FontSize',12,'FontWeight','bold')
group_legend={'MPAS-O hex ','POP z-level','POP z-star',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5 3.5])

fig=[char(dir) '_1km_dRPEdt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 4 km 

dir={'m77'}; abc='fghijFGHIJ'; 

nu_h=[1 5 10 20 200];
grid_spacing=4e3;

hist_files = {...
''};

time_fields=50*[0:6]; % time in days
min_n=1;max_n = length(time_fields);

kmt_file = 'kmt_BCC_4kms.da';
grid_file = 'grid_BCC_4kms.da';
depth_file = 'in_depths.dat';

title_txt={
    'baroclinic channel',...
	  };
if (1==1)
rpeNorm = zeros(length(time_fields),length(abc));
rpeTot = zeros(length(time_fields),length(abc));
keMeanVolume = zeros(length(time_fields),length(abc));
meanDrpeDt4kmPOP100day = zeros(1,length(abc));
keMeanTime4kmPOP = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
[time,rpeTot(:,j),rpeNorm(:,j),DrpeDt,meanDrpeDt4kmPOP(j),keMeanTime4kmPOP(j),keMeanVolume(:,j)] ...
    = sub_rpe(wd,dir,abc(j),hist_files, kmt_file, grid_file, depth_file, ...
    time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt4km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt4kmPOP100day(j));
end

end

%%%%%% plot MPAS dRPE/dt versus grid Reynolds number

figure(1); clf;
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime4km);
for j=2:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt4km320day(iBeg:iEnd),'k');
  set(h,'Marker',marker(j),'LineWidth',1)
  if j~=2; set(h,'LineStyle','none'); end
  grid on
  %axis([1e-1 2e3 1e-6 9e-3])
  hold on
end

%%%%%% plot MPAS dRPE/dt versus grid Reynolds number QUADS

load('/home/mpeterse/mpas/matlab/rpe/m83_bcl_ch_quads/data/m83QuadHex4km.mat')

m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTimeQuadHex4km);
for j=1:1 %length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDtQuadHex4km(iBeg:iEnd),'k');
  set(h,'Color','k','MarkerFaceColor','k');
  set(h,'Marker','s','LineWidth',1)
  set(h,'LineStyle','none')
  hold on
end

%%%%%% plot POP dRPE/dt versus grid Reynolds number

figure(1); %clf;
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime4kmPOP);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt4kmPOP(iBeg:iEnd),'r');
  set(h,'Marker',marker(j),'LineWidth',1)
  if j~=2; set(h,'LineStyle','none'); end
  grid on
  %axis([1e-1 2e3 1e-6 9e-3])
  hold on
end

h=loglog(mitgcm_4km_Re,mitgcm_4km_drpe_dt,'^b');
set(h,'MarkerFaceColor','g','Color','g');
h=loglog(mom_4km_Re,mom_4km_drpe_dt,'om');
set(h,'MarkerFaceColor','b','Color','b');

%title('MPAS-O Baroclinic channel, 4km')
grid on
axis([.7 1e3 1e-6 3e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
h=text(1,1.5e-3,'4 km resolution')
set(h,'FontSize',12,'FontWeight','bold')
group_legend={'MPAS-O hex','MPAS-O quad','POP z-level','POP z-star',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])

fig=[char(dir) '_4km_dRPEdt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 10 km 

dir={'m77'}; abc='abcdeABCDE'; 

nu_h=[1 5 10 20 200];
grid_spacing=10e3;

hist_files = {...
''};

time_fields=50*[0:6]; % time in days
min_n=1;max_n = length(time_fields);

kmt_file = 'kmt_BCC_10km.da';
grid_file = 'grid_BCC_10km.da';
depth_file = 'in_depths.dat';

title_txt={
    'baroclinic channel',...
	  };
if (1==1)
rpeNorm = zeros(length(time_fields),length(abc));
rpeTot = zeros(length(time_fields),length(abc));
keMeanVolume = zeros(length(time_fields),length(abc));
meanDrpeDt10kmPOP100day = zeros(1,length(abc));
keMeanTime10kmPOP = zeros(1,length(abc));
vertTransportMean = zeros(1,length(abc));
nVertLevels = 20;
vertTransportZ = zeros(nVertLevels,length(abc));
for j=1:length(abc)
[time,rpeTot(:,j),rpeNorm(:,j),DrpeDt,meanDrpeDt10kmPOP(j),keMeanTime10kmPOP(j),keMeanVolume(:,j)] ...
    = sub_rpe(wd,dir,abc(j),hist_files, kmt_file, grid_file, depth_file, ...
    time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt10km100day ' char(dir) abc(j) ': %e \n'],meanDrpeDt10kmPOP100day(j));
end

end

%%%%%% plot MPAS dRPE/dt versus grid Reynolds number

figure(1); clf;
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime10km);
for j=2:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt10km320day(iBeg:iEnd),'k');
  set(h,'Marker',marker(j),'LineWidth',1)
  if j~=2; set(h,'LineStyle','none'); end
  grid on
  %axis([1e-1 2e3 1e-6 9e-3])
  hold on
end

%%%%%% plot MPAS dRPE/dt versus grid Reynolds number QUADS

load('/home/mpeterse/mpas/matlab/rpe/m83_bcl_ch_quads/data/m83QuadHex10km.mat')

m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTimeQuadHex10km);
for j=1:1 %length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDtQuadHex10km(iBeg:iEnd),'k');
  set(h,'Color','k','MarkerFaceColor','k');
  set(h,'Marker','s','LineWidth',1)
  set(h,'LineStyle','none')
  hold on
end

%%%%%% plot POP dRPE/dt versus grid Reynolds number

figure(1); %clf;
m=length(nu_h); % # experiments in set
vel_scale = sqrt(2*keMeanTime10kmPOP);
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,meanDrpeDt10kmPOP(iBeg:iEnd),'r');
  set(h,'Marker',marker(j),'LineWidth',1)
  if j~=2; set(h,'LineStyle','none'); end
  grid on
  %axis([1e-1 2e3 1e-6 9e-3])
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
h=text(1,1.5e-3,'10 km resolution')
set(h,'FontSize',12,'FontWeight','bold')
group_legend={'MPAS-O hex','MPAS-O quad','POP z-level','POP z-star',...
	      'MITGCM','MOM',
	     };
legend(group_legend,'location','SouthEast')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])

fig=[char(dir) '_10km_dRPEdt'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return
%%%%%% plot RPE norm versus time

figure(3); clf;
m=length(nu_h); % # experiments in set
for j=1:length(abc)
  h=plot(time_fields,rpeTot(:,j));
  set(h,'LineWidth',1,'Color',cl(j,:))
  if j>m; set(h,'LineStyle','--'); end
  grid on
  %axis([0 300 0  14e-7])
  hold on
end
title('POP z-level, POP z-star (dashed) baroclinic channel, 10km')
xlabel('time, days')
ylabel('RPE, J/m^2')

group_legend={'\nu_h=1','\nu_h=5','\nu_h=10','\nu_h=20','\nu_h=200'...
	      };
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])

fig=[char(dir) '_10km_RPEtot'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return



%%%%%% plot RPE norm versus time
 
figure(12); clf;
m=length(nu_h); % # experiments in set
for j=1:length(abc)
  h=plot(time_fields,rpeNorm(:,j));
  set(h,'LineWidth',1,'Color',cl(j,:))
  if j>m; set(h,'LineStyle','--'); end
  grid on
  axis([0 300 0  14e-7])
  hold on
end
title('POP z-level, POP z-star (dashed) baroclinic channel, 1km')
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')

group_legend={'\nu_h=1','\nu_h=5','\nu_h=10','\nu_h=20','\nu_h=200'...
	      };
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])
fig=[char(dir) '_1km_RPEnorm'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);


%%%%%% plot RPE norm versus time

figure(13); clf;
m=length(nu_h); % # experiments in set
for j=1:length(abc)
  h=plot(time_fields,rpeTot(:,j));
  set(h,'LineWidth',1,'Color',cl(j,:))
  if j>m; set(h,'LineStyle','--'); end
  grid on
  %axis([0 300 0  14e-7])
  hold on
end
title('POP z-level, POP z-star (dashed) baroclinic channel, 1km')
xlabel('time, days')
ylabel('RPE, J/m^2')

group_legend={'\nu_h=1','\nu_h=5','\nu_h=10','\nu_h=20','\nu_h=200'...
	      };
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])

fig=[char(dir) '_1km_RPEtot'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return

%%%%%% plot RPE norm versus time

figure(2); clf;
m=length(nu_h); % # experiments in set
for j=1:length(abc)
  h=plot(time_fields,rpeNorm(:,j));
  set(h,'LineWidth',1,'Color',cl(j,:))
  if j>m; set(h,'LineStyle','--'); end
  grid on
  axis([0 300 0  14e-7])
  hold on
end
title('POP z-level, POP z-star (dashed) baroclinic channel, 4km')
xlabel('time, days')
ylabel('(RPE-RPE(0))/RPE(0)')

group_legend={'\nu_h=1','\nu_h=5','\nu_h=10','\nu_h=20','\nu_h=200'...
	      };
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])
fig=[char(dir) '_4km_RPEnorm'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);


%%%%%% plot RPE norm versus time

figure(3); clf;
m=length(nu_h); % # experiments in set
for j=1:length(abc)
  h=plot(time_fields,rpeTot(:,j));
  set(h,'LineWidth',1,'Color',cl(j,:))
  if j>m; set(h,'LineStyle','--'); end
  grid on
  %axis([0 300 0  14e-7])
  hold on
end
title('POP z-level, POP z-star (dashed) baroclinic channel, 4km')
xlabel('time, days')
ylabel('RPE, J/m^2')

group_legend={'\nu_h=1','\nu_h=5','\nu_h=10','\nu_h=20','\nu_h=200'...
	      };
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 5.0 3.5])

fig=[char(dir) '_4km_RPEtot'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

return


