% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:33];
dir={'m62'}; abc='ghijk'; dims=1*[40 143]; % change after m51k is done
nu_h=[1 5 10 20 200];
grid_spacing=4e3;

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'baroclinic channel',...
	  };

min_n=2;
max_n = 33;

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

%%%%%%%%%%%%%%%%%%%%%
tau_Dlf = [.1 1 5 10 100 1000];

vel_scale = sqrt(2*keMeanTime);
% for dam break, velocity scale is from eqn 5:
% or:
%rho1=999;
%rho2=997;
%vel_scale = 1/2*sqrt(9.8*500*(rho1-rho2)/((rho1+rho2)/2)) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, many viscosities

group_legend={'\nu_h=1','\nu_h=5','\nu_h=10','\nu_h=20','\nu_h=200','\nu_h=800',...
	      };
figure(11); clf
nsims=5;
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
%title('MPAS-O Baroclinic channel, 10km')
axis tight
legend(group_legend,'location','NorthWest')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4 3])

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
%%%%%% plot dRPE/dt versus grid Reynolds number

figure(20); clf
m=length(nu_h); % # experiments in set
marker='s*';
for j=length(abc)/m:-1:1
  iBeg=(j-1)*m+1
  iEnd=j*m
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h
  h=loglog(gridRe,meanDrpeDt(iBeg:iEnd),'-k');
  %set(h,'Color',cl(j,:))
  set(h,'Marker',marker(j))
  hold on
end
grid on
axis([.1 1e4 9e-6 3e-3])

load '~/mpas/email/130911_Ilicak_fig_data/fig17.mat'
h=loglog(mitgcm_10km_Re,mitgcm_10km_drpe_dt,'^k');
set(h,'MarkerFaceColor','k');
%set(h,'Color',   [ 1    0    0],'LineWidth',1) ;% red
h=loglog(mom_10km_Re,mom_10km_drpe_dt,'ok');
set(h,'MarkerFaceColor','k');
%set(h,'Color',   [ 0    0.6  0.9],'LineWidth',1) ;% blue
%j=j+1; set(h,'Color',cl(j,:))

%title('MPAS-O Baroclinic channel, 10km')
xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
group_legend={'MPAS-O, z-level','MPAS-O, z-star','MITGCM','MOM',
	     };
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


