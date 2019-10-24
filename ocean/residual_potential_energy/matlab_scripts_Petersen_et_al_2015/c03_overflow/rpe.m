% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:21];
dir={'c03'}; abc='bcdefg'; dims=1*[4 230];
dir={'c03'}; abc='hijklm'; dims=1*[4 230];

time_fields=[1:41];
dir={'c05'}; abc='no'; dims=1*[4 230];
min_n = 14;
max_n = 38;

time_fields=[1:8];
dir={'c05'}; abc='jklmno'; dims=1*[4 230];
min_n = 1;
max_n = 7;

time_fields=[1:40];
dir={'c18'}; abc='fglm'; dims=1*[4 230];
min_n = 14;
max_n = 38;

time_fields=[1:13];
dir={'c18'}; abc='cdefgijklmopqrs'; dims=1*[4 230];
min_n = 1;
max_n = 13;

time_fields=[1:38];
dir={'c18'}; abc='nopqrs'; dims=1*[4 230];
min_n = 14;
max_n = 38;

time_fields=[1:38];
dir={'c19'}; abc='hijklmfg'; dims=1*[4 230];
nu_h=[10.^[-2:3] 10.^[2:3] 10.^[-2:3]];
min_n = 14;
max_n = 38;

time_fields=[1:38];
dir={'c18'}; abc='hijklmfgnopqrs'; dims=1*[4 230];
nu_h=[10.^[-2:3] 10.^[2:3] 10.^[-2:3]];
min_n = 14;
max_n = 38;

%group_legend={'MPAS-O, z-full cell','MPAS-O, z-pbc','MPAS-O, \sigma','MITGCM','MOM'};
%group_legend={'MPAS-O, z*-full cell','MPAS-O, z*-pbc','MPAS-O, \sigma','MPAS-O, z-full cell','MPAS-O, z-pbc'};
group_legend={'MPAS-O, z*-full cell','MPAS-O, z*-pbc','MPAS-O, \sigma'}


netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'overflow',...
	  };

%load('data/m66_overflow_131021.mat')
if (1==2)
load '~/mpas/email/130911_Ilicak_fig_data/fig10.mat'

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
save(['data/' char(dir) abc '_dRPEdt.mat'])
end

%load(['data/' char(dir) abc '_dRPEdt.mat'])


cl=[...
  0    0    0   ; ... % black
  0    0.6  0.9 ; ... % blue
  1    0    0   ; ... % red
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

%%%%%%% simple plot of rpe

figure(20)
h=plot(time/3600,rpeNorm(:,:),'-');
grid on
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')
%axis([0 10 -1e-5 12e-5])
group_legend={'nu_h=0.01','nu_h=0.1','nu_h=1','nu_h=10','nu_h=100','nu_h=1000'};
%group_legend={'nu_h=100','nu_h=1000'};
legend(group_legend,...
       'location','NorthWest')
axis tight
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])

fig=[char(dir) abc '_rpe_short'];
print('-djpeg',['f/' fig '.jpg']);
print('-depsc2',['f/' fig '.eps']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

%%%%%%%%%%%%%%%%%%%%%
%nu_h=[10.^[-2:3]];
tau_Dlf = [.1 1 5 10 100 1000];
grid_spacing=1000;

%vel_scale = sqrt(2*keMeanTime);
% for overflow, velocity scale is from eqn 5:
% or:
rho1=999;
rho2=997;
vel_scale = 1/2*sqrt(9.8*500*(rho1-rho2)/((rho1+rho2)/2)) 



%%%%%% plot dRPE/dt versus grid Reynolds number

figure(4); clf
m=length(nu_h); % # experiments in set
marker='s*^';
iBeg=[1 7 9]
iEnd=[6 8 14]
for j=1:3 %length(abc)/m
  %iBeg=(j-1)*m+1;
  %iEnd=j*m;
  gridRe = grid_spacing*vel_scale./nu_h(iBeg(j):iEnd(j))
  h=loglog(gridRe,meanDrpeDt(iBeg(j):iEnd(j)),'-k');
  set(h,'Color',cl(j,:),'Marker',marker(j))
  hold on
end
grid on

h=loglog(mitgcm_Re,mitgcm_drpe_dt,'^g');
set(h,'MarkerFaceColor','g','MarkerSize',10);
h=loglog(mom_Re,mom_drpe_dt,'ob');
set(h,'MarkerFaceColor','b');

% put partial cell marker on top
for j=2 %length(abc)/m
  %iBeg=(j-1)*m+1;
  %iEnd=j*m;
  gridRe = grid_spacing*vel_scale./nu_h(iBeg(j):iEnd(j))
  h=loglog(gridRe,meanDrpeDt(iBeg(j):iEnd(j)),'-k');
  set(h,'Color',cl(j,:),'Marker',marker(j))
end

xlabel('grid Reynolds number')
ylabel('dRPE/dt, W/m^2')
group_legend={'MPAS-O z* full cell','MPAS-O z* partial cell','MPAS-O sigma',...
	      'MITgcm','MOM'}
legend(group_legend,...
       'location','NorthWest')
axis([1 1e6 .08 5])
set(gca,'XTick',10.^[0:6])
set(gca,'XMinorGrid','off')

set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4 3])
fig=[char(dir) abc '_DrpeDt_t37'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
%fig=[char(dir) abc '_DrpeDt'];

return
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
