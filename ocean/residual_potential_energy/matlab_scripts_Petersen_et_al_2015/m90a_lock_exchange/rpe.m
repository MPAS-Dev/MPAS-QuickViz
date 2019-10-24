% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:121];
dir={'m46'}; abc='klmnopabcdef'; dims=1*[2 145];

time_fields=[1:18];
dir={'m90'}; abc='a'; dims=1*[4 148];

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'dam break',...
	  };

min_n=17; % time 16:50
max_n = 18; % time 17:10

if (1==1)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, one only

group_legend={'MPAS-Ocean, \Delta z=1'
	      };
figure(12); clf
nsims=1;
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot(time/3600,rpeNorm(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
end

load '~/mpas/email/130911_Ilicak_fig_data/fig4.mat'
h=loglog(mitgcm_delta_x_500_time_days,mitgcm_delta_x_500_drpe_dt,'-');
set(h,'Color',   [ 1    0    0],'LineWidth',1) ;% red
h=loglog(mom_delta_x_500_time_days,mom_delta_x_500_drpe_dt,'-');
set(h,'Color',   [ 0    0.6  0.9],'LineWidth',1) ;% blue
h=loglog(roms_delta_x_500_time_days,roms_delta_x_500_drpe_dt,'-');
set(h,'Color',   [ 0.1  1    0.1],'LineWidth',1) ;% green
h=loglog(gold_delta_z_10_time_days,gold_delta_z_10_drpe_dt,'-');
set(h,'Color',[255, 215, 0]/255,'LineWidth',1); % dark yellow
%j=j+1; set(h,'Color',cl(j,:))


grid on
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')
axis([0 17 -1e-5 8e-5])
%axis tight
group_legend={...
 'MPAS-O','MITGCM','MOM','ROMS','GOLD'};
%title('lock exchange, \nu_h=0.01')
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4 3])

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

%%%%%%%%%%%%%%%%%%%%%
group_legend={'MPAS-O: z-level','MPAS-O: z-star','MITGCM'... 
	      };
nu_h=[10.^[-2:2] 200];
tau_Dlf = [.1 1 5 10 100 1000];
grid_spacing=500;

%vel_scale = sqrt(2*keMeanTime);
% for dam break, velocity scale is from eqn 5:
% vel_scale = (50e3-32e3)/10/3600; % from fig 3b, gives 0.5 m/s
% or:
rho1=1000;
rho2=995;
vel_scale = 1/2*sqrt(9.8*20*(rho1-rho2)/((rho1+rho2)/2)) ;
% previous gives 0.4956 m/s, same as from fig 3b!

%%%%%% plot dRPE/dt versus grid Reynolds number

figure(10); clf
m=length(nu_h); % # experiments in set
marker='s*';
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale./nu_h;
  h=loglog(gridRe,meanDrpeDt(iBeg:iEnd),'-*k');
  %set(h,'Color',cl(j,:))
  set(h,'Marker',marker(j))
  set(h,'LineWidth',1)
  hold on
end

load '~/mpas/email/130911_Ilicak_fig_data/fig6.mat'
h=loglog(mitgcm_dx_500_nu_changes_Re,mitgcm_dx_500_nu_changes_drpe_dt ,'^r-');
  set(h,'LineWidth',1)
%j=j+1; set(h,'Color',cl(j,:))

grid on
axis([1 1e5 1e-4 5e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
%title('lock exchange')
legend(group_legend,'location','SouthEast')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4 3])
fig=[char(dir) '_DrpeDt'];
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

figure(11); clf
nsims=6;
for j=1:nsims
  %h=plot(time/3600,slidefilter(rpeNorm(:,j),8),'-');
  h=plot(time/3600,rpeNorm(:,j),'-');
  set(h,'Color',cl(j,:),'LineWidth',1)
  hold on
end
grid on
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')
%title('lock exchange')
axis([0 17 -1e-5 8e-5])
%axis tight
group_legend={'\nu_h=0.01','\nu_h=0.1','\nu_h=1','\nu_h=10','\nu_h=100','\nu_h=200',...
	      };
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% plot vert transport versus grid Reynolds number

figure(13); clf
m=length(nu_h); % # experiments in set
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,vertTransportZ(nVertLevels/2,iBeg:iEnd),'-*');
  set(h,'Color',cl(j,:))
  hold on
end
grid on
%axis([.1 1e5 5e-5 1e-3])
xlabel('grid Reynolds number')
ylabel('mean(abs(vert transport)) at mid-depth m/s')
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[1 3 4 3])
fig=[char(dir) '_vertTrans'];
print('-depsc2',['f/' fig '.eps']);
unix(['convert f/' fig '.eps f/' fig '.pdf']);
subplot('position',[0 .95 1 .05]); axis off
h=text(.5,.4,title_txt);
set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
text(.005,.7,[ date ]);
text(.005,0,[char(dir)]);% abc]);
%fig=[char(dir) abc '_DrpeDt'];

%%%%%% plot vert transport versus grid Reynolds number

figure(14); clf

m=length(nu_h); % # experiments in set
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  h=semilogx(vertTransportZ(:,iEnd),1:nVertLevels,'-*');
  set(h,'Color',cl(j,:))
  hold on
end
set(gca,'YDir','reverse')
grid on
axis([1e-6 4e-4 1 20])
ylabel('vertical layer')
xlabel('mean(abs(vert transport)), m/s')
legend(group_legend,'location','SouthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[1 3 6.0 4.5])
%fig=[char(dir) abc '_DrpeDt'];
fig=[char(dir) '_vertTransZ'];
print('-depsc2',['f/' fig '.eps']);
unix(['convert f/' fig '.eps f/' fig '.pdf']);
subplot('position',[0 .95 1 .05]); axis off
h=text(.5,.4,title_txt);
set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
text(.005,.7,[ date ]);
text(.005,0,[char(dir)]);% abc]);



return
figure(11); clf
group_legend={'\nu_h=0.01','\nu_h=0.1','\nu_h=1','\nu_h=10','\nu_h=100',...
	      };
for j=1:4
  h=plot(time/3600,slidefilter(rpeNorm(:,j),12));
  set(h,'Color',cl(j,:))
  hold on
end
grid on
xlabel('time, hours')
ylabel('(rpe-rpe(0))/rpe(0)')
%axis([0 200 0 5e-6])
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[1 3 6.0 4.5])
subplot('position',[0 .95 1 .05]); axis off
h=text(.5,.4,title_txt);
set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
text(.005,.7,[ date ]);
text(.005,0,[char(dir) abc]);
fig=[char(dir) abc '_rpe'];
print('-depsc2',['f/' fig '.eps']);
unix(['convert f/' fig '.eps f/' fig '.pdf']);



% drpe/dt
nu_h=10.^[-2:2];
grid_spacing=5e3;
vel_scale = sqrt(2*keMeanTime);
gridRe = grid_spacing*vel_scale./nu_h;

figure(12); clf
loglog(gridRe,meanDrpeDt,'-*')
grid on
axis([1 1e5 1e-4 1e-3])
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
legend('z-star',...
       'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[1 3 6.0 4.5])
subplot('position',[0 .95 1 .05]); axis off
h=text(.5,.4,title_txt);
set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
text(.005,.7,[ date ]);
text(.005,0,[char(dir) abc]);
fig=[char(dir) abc '_DrpeDt'];
print('-depsc2',['f/' fig '.eps']);
unix(['convert f/' fig '.eps f/' fig '.pdf']);
