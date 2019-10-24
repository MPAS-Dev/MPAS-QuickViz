% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:201];
%dir={'m52'}; abc='abcdefghijklmnopqrstuvwx'; dims=1*[2 57];
dir={'m52'}; abc='abcdefghijklmnopqrstuvwx'; dims=1*[2 57];
grid_spacing=5e3;
nu_h=[.01 1 15 150];
% load('data/m52_igw_m52a-x_vert_transport.mat')

netcdf_file = 'output.0000-01-01_00.00.00.nc';

title_txt={
    'internal wave test',...
	  };

min_n=11;
max_n = 101;

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

vel_scale = sqrt(2*keMeanTime);

cl=[...
  0    0    0   ; ... % black
  0    0    0   ; ... % black
  0    0.6  0.9 ; ... % blue
  1    0    0   ; ... % red
  0    0.4  0.1 ; ... % dark green
  0.1  1    0.1 ; ... % green
];
mk = 's**********';
blue = [0    0.6  0.9];
%%%%%%%%%%%%%%%%%%%%%
group_legend={'z-level','z-star',... 
	       '\tau_{Dlf}=0.1',...
	       '\tau_{Dlf}=1',...
	       '\tau_{Dlf}=10',...
	       '\tau_{Dlf}=100',...
	       '\tau_{Dlf}=1000',...
	      };
tau_Dlf = [.1 1 5 10 100 1000];

% for dam break, velocity scale is from eqn 5:
% or:
%rho1=999;
%rho2=997;
%vel_scale = 1/2*sqrt(9.8*500*(rho1-rho2)/((rho1+rho2)/2)) 

%%%%%% plot vert transport versus z
figure(14); clf

m=length(nu_h); % # experiments in set
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  h=semilogx(vertTransportZ(:,iBeg),[0:25:499],'-*');
  set(h,'Color',cl(j,:),'Marker',mk(j),'LineWidth',1)
  hold on
end
set(gca,'YDir','reverse')
grid on
axis([.99e-5 3e-3 0 500])
ylabel('depth, m')
xlabel('rms of vertical transport, m/s')
legend(group_legend,'location','West')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 6 4])
%fig=[char(dir) abc '_DrpeDt'];
fig=[char(dir) '_vertTransZ'];
print('-depsc2',['f/' fig '.eps']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
%subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir)]);% abc]);


%%%%%% plot vert transport versus grid Reynolds number

figure(13); clf
m=length(nu_h); % # experiments in set
for j=1:length(abc)/m
  iBeg=(j-1)*m+1;
  iEnd=j*m;
  gridRe = grid_spacing*vel_scale(iBeg:iEnd)./nu_h;
  h=loglog(gridRe,vertTransportZ(nVertLevels/2,iBeg:iEnd),'-');
  set(h,'Color',cl(j,:),'Marker',mk(j),'LineWidth',1)
  hold on
end
grid on
axis([1 1e5 9.9e-6 4e-3])
xlabel('grid Reynolds number')
ylabel('rms of vertical transport at mid-depth, m/s')
legend(group_legend,'location','EastOutside')
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[1 3 7 4])
fig=[char(dir) '_vertTrans'];
print('-depsc2',['f/' fig '.eps']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);
%subplot('position',[0 .95 1 .05]); axis off
%h=text(.5,.4,title_txt);
%set(h,'HorizontalAlignment','center','FontWeight','bold','FontSize',14)
%text(.005,.7,[ date ]);
%text(.005,0,[char(dir)]);% abc]);
%fig=[char(dir) abc '_DrpeDt'];
