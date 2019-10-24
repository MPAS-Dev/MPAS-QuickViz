function [time,rpe,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime,vertTransportMean,vertTransportMeanZ, ...
	headLocation,DrpeDx ] ...
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

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xCell'));
nCells = size(work,1);
xCell = reshape(work(:,1), dims);

yCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yCell'));
yCell = reshape(yCellFull(:,1), dims);

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yEdge'));
yEdge = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xEdge'));
xEdge = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'areaCell'));
areaCell = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'maxLevelCell'));
maxLevelCell = squeeze(work(:,1));

refLayerThickness = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refLayerThickness'));
refBottomDepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
refZBot = -refBottomDepth;

hFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'layerThickness'));
densityFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'density'));
kineticEnergyCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kineticEnergyCell'));
vertTransportFull = abs(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'vertTransportVelocityTop')));
kmax = size(hFull,1);

bottomMax = max(refBottomDepth);
gravity = 9.80616;

time = zeros(size(time_fields));
rpe = zeros(size(time_fields));
headLocation = zeros(size(time_fields));
keMeanVolume = zeros(size(time_fields));
vertTransportVolume = zeros(size(time_fields));
vertTransportVolumeZ = zeros(kmax,length(time_fields));

% note: time_fields(1) needs to be IC time to compute volume correctly.
  h = squeeze(hFull(:,:,time_fields(1)));  
  layerVolume = zeros(kmax,1);
  for iCell=1:nCells
    for k=1:maxLevelCell(iCell)
      layerVolume(k) = layerVolume(k) + h(k,iCell)*areaCell(iCell);
    end
  end 
layerArea = layerVolume./refLayerThickness;

fprintf(['begin sort, ' char(dir) abc ' nt='])
for nt=1:length(time_fields)
  fprintf([' ' num2str(nt)])
  t = xtime(:,time_fields(nt))';  % time string
  time(nt)=str2num(t(18:19)) ... % s
         + str2num(t(15:16))   *60 ... %  min
         + str2num(t(12:13))   *60*60 ... %  hour
         +(str2num(t( 9:10))-1)*60*60*24 ... %  day
         +(str2num(t( 6: 7))-1)*60*60*24*30 ... %  month
         + str2num(t( 1: 4))   *60*60*24*30*360;  %  year

  h = squeeze(hFull(:,:,time_fields(nt)));  
  density = squeeze(densityFull(:,:,time_fields(nt)));

  headLocationZ = zeros(1,kmax);
  headDensity=998;
  for k=1:kmax
    temp = squeeze(densityFull(k,:,time_fields(nt)));
    if max(temp)>headDensity
      denseWaterIndices = find(temp>headDensity);
      headLocationZ(k) = max(yCellFull(denseWaterIndices,1));
    end
  end
  headLocation(nt) = max(headLocationZ);
  
  vol_1D = zeros(1,kmax*nCells);
  density_1D = zeros(1,kmax*nCells);
  j=0;
  for iCell=1:nCells
    for k=1:maxLevelCell(iCell)
      j = j+1;
      vol_1D(j) = h(k,iCell)*areaCell(iCell);
      density_1D(j) = density(k,iCell);
    end
  end 
  nCells_1D = j;

  [sortedDensity,sortedIndex] = sort( density_1D(1:nCells_1D) );
  sortedVolume = vol_1D(sortedIndex);
  %fprintf('  done \n')

  sortedRPE = zeros(1,nCells_1D);
  sortedZMid = zeros(1,nCells_1D);
  zTop = 0.0;
  k=1; % start at top layer
  for j=1:nCells_1D
    deltaZ = sortedVolume(j)/layerArea(k);
    % try taking this out.  That means topography is never
    % considered, use surface area all the way down.  layerArea(k)
    % is always k=1, so surface.  No, didn't make any difference
    if zTop-deltaZ<refZBot(k)&k<kmax
      upperVolume = layerArea(k)*(zTop - refZBot(k));
      lowerVolume = sortedVolume(j) - upperVolume;
      deltaZLower = lowerVolume/layerArea(k+1);
      zBot = refZBot(k) - deltaZLower;
      k=k+1;
    else
      zBot = zTop - deltaZ;
    end
    sortedZMid(j) = 0.5*(zTop+zBot);
    sortedRPE(j) = gravity*sortedDensity(j)*(sortedZMid(j)+bottomMax)*sortedVolume(j);
    zTop = zBot;
  end

  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpe(nt) = sum(sortedRPE)/sum(areaCell);

  % note: this only works with no land cells
  keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_fields(nt)))));
  vertTransportVolume(nt) = mean(mean((vertTransportFull(:,:, ...
						  time_fields(nt)))));
  for k=1:kmax
    vertTransportVolumeZ(k,nt) = mean(mean((vertTransportFull(k,:,time_fields(nt)))));
  end
end
fprintf('\n')
rpe(1);

rpeNorm = (rpe - rpe(1))/rpe(1);

DrpeDt = derivative_ord2(time,rpe);
DrpeDx = derivative_ord2(headLocation,rpe);

n=length(DrpeDt);

meanDrpeDt = mean(DrpeDt(min_n:max_n));
keMeanTime = mean(keMeanVolume(min_n:max_n));

vertTransportMean = mean(vertTransportVolume(min_n:max_n));
  for k=1:kmax
    vertTransportMeanZ(k) = mean((vertTransportVolumeZ(k,:)));
  end
  
