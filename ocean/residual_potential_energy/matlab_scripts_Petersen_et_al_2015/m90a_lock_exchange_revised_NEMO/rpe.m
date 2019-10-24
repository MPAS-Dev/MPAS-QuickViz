% script rpe
% compute resting potential energy, ala Ilicak ea 2012
% Mark Petersen, LANL, Jan 2013

% working directory, where data is kept:
%wd = '/local1/mpetersen/runs/';
wd = '/var/tmp/mpeterse/runs/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_fields=[1:18];
dir={'m90'}; abc='a'; dims=1*[4 148];
netcdf_file = 'output.0000-01-01_00.00.00.nc';
filename = [wd char(dir) abc '/' netcdf_file];

min_n=17; % time 16:50
max_n = 18; % time 17:10

[time,rpeTot,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime] ...
    = sub_rpe(filename,...
    dims,time_fields,min_n,max_n,char(title_txt));
  fprintf(['meanDrpeDt ' char(dir) abc ': %e \n'],meanDrpeDt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RPE with time, one only

group_legend={'MPAS-Ocean, \Delta z=1'
	      };
figure(12); clf

h=plot(time/3600,rpeNorm,'-');

grid on
xlabel('time, hours')
ylabel('(RPE-RPE(0))/RPE(0)')
axis([0 17 -1e-5 8e-5])
