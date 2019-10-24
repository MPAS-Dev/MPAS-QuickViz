function [time,rpeTot,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime,vertTransportMean,vertTransportMeanZ] ...
    = sub_rpe( ...
    wd,dir,abc,netcdf_file, ...
    dims,time_fields,min_n,max_n,title_txt)

% Compute Resting Potential Energy
% Mark Petersen, LANL, Jan 2013

px = [.53];
py=[.52]; % Midpoint position of plots
pw = [.87];  ph=[.75]; % width and height of plots

filename = [wd char(dir) abc '/' netcdf_file];
ncid = netcdf.open(filename,'nc_nowrite');

xtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xtime'));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'areaCell'));
areaCell = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'maxLevelCell'));
maxLevelCell = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
refBottomDepth = squeeze(work(:,1));

hFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'layerThickness'));
densityFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'potentialDensity'));
temperatureFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'temperature'));
kineticEnergyCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kineticEnergyCell'));
%vertTransportFull = abs(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'vertTransportVelocityTop')));
K = size(hFull,1);
nCells = size(areaCell,1);
maxBottom = max(refBottomDepth);
gravity = 9.80616;

time = zeros(size(time_fields));
rpe = zeros(size(time_fields));
keMeanVolume = zeros(size(time_fields));
%vertTransportVolume = zeros(size(time_fields));
%vertTransportVolumeZ = zeros(K,length(time_fields));

for nt=1:length(time_fields)
  t = xtime(:,time_fields(nt))';  % time string
  time(nt)=str2num(t(18:19)) ... % s
         + str2num(t(15:16))   *60 ... %  min
         + str2num(t(12:13))   *60*60 ... %  hour
         +(str2num(t( 9:10))-1)*60*60*24 ... %  day
         +(str2num(t( 6: 7))-1)*60*60*24*30 ... %  month
         + str2num(t( 1: 4))   *60*60*24*30*360;  %  year

  h = squeeze(hFull(:,:,time_fields(nt)));  
  density = squeeze(densityFull(:,:,time_fields(nt)));

%%%%%%%%%%%%%% test  
%    if (linear_eos)
  TEMP = squeeze(temperatureFull(:,:,time_fields(nt)));
        config_eos_linear_alpha = 0.2;
        config_eos_linear_beta = 0.8;
        config_eos_linear_Tref = 5.0;
        config_eos_linear_Sref = 35.0;
        config_eos_linear_densityref = 1000.0;
     density = config_eos_linear_densityref - config_eos_linear_alpha*(TEMP-config_eos_linear_Tref);
%  end
%%%%%%%%%%%%%% test end
  
  vol_1D = zeros(1,K*nCells);
  density_1D = zeros(1,K*nCells);
  i=0;
  for iCell=1:nCells
    for k=1:maxLevelCell(iCell)
      i = i+1;
      vol_1D(i) = h(k,iCell)*areaCell(iCell);
      density_1D(i) = density(k,iCell);
    end
  end 
  nCells_1D = i;

%  fprintf('sorting...')
  [density_sorted,sorted_ind] = sort( density_1D(1:nCells_1D) );
%  fprintf('done \n')
  vol_sorted = vol_1D(sorted_ind);

  zMid = zeros(1,nCells_1D);
  rpeSorted = zeros(1,nCells_1D);
  z = maxBottom;
  % this only works for flat-bottom cases:
  sumCellArea = sum(areaCell);
  for i=1:nCells_1D
    thickness = vol_sorted(i)/sumCellArea;
    zMid(i) = z - thickness/2;
    z = z - thickness;
    rpeSorted(i) = gravity*density_sorted(i)*zMid(i)*vol_sorted(i);
  end
  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpeTot(nt) = sum(rpeSorted)/sum(areaCell);

  % note: this only works with no land cells
  keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_fields(nt)))));
%  vertTransportVolume(nt) = mean(mean((vertTransportFull(:,:, ...
%						  time_fields(nt)))));
%  for k=1:K
%    vertTransportVolumeZ(k,nt) = mean(mean((vertTransportFull(k,:,time_fields(nt)))));
%  end
end
rpeTot(1);

rpeNorm = (rpeTot - rpeTot(1))/rpeTot(1);

DrpeDt = derivative_ord2(time,rpeTot);

n=length(DrpeDt);
%%%%%%%%%%% make sure this is right
%min_n=1;
%min_n=floor(n/2);
%max_n = floor(n/2);
%max_n = n;

%min_n=90;
%max_n=96;

meanDrpeDt = mean(DrpeDt(min_n:max_n));
keMeanTime = mean(keMeanVolume(min_n:max_n));

%vertTransportMean = mean(vertTransportVolume(min_n:max_n));
%  for k=1:K
%    vertTransportMeanZ(k) = mean((vertTransportVolumeZ(k,:)));
%  end

vertTransportMean=0;
vertTransportMeanZ=0;
