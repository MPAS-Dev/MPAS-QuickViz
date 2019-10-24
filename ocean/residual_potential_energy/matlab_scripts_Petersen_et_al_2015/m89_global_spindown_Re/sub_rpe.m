function [time,rpeTot,rpeNorm,DrpeDt,keMeanVolume, ... %vertTransportVolume,vertTransportMeanZ, ...
	  rpeMean,rpePert] ...
    = sub_rpe( ...
    wd,dir,abc,netcdf_file, ...
    time_fields,min_n,max_n,title_txt,multiple_files)

% Compute Resting Potential Energy
% Mark Petersen, LANL, Jan 2013

if (multiple_files==true)
  filename = [wd char(dir) abc '/' netcdf_file '_time_' num2str(time_fields(1)) '.nc'];
  ncid = netcdf.open(filename,'nc_nowrite');
else
  filename = [wd char(dir) abc '/' netcdf_file ]
  ncid = netcdf.open(filename,'nc_nowrite');
  xtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xtime'));
  hFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'layerThickness'));
  densityFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'potentialDensity'));
  kineticEnergyCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kineticEnergyCell'));
end

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'areaCell'));
areaCell = squeeze(work(:,1));
nCells = length(areaCell);

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'maxLevelCell'));
maxLevelCell = squeeze(work(:,1));
refLayerThickness = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refLayerThickness'));
nVertLevels = size(refLayerThickness,1);
refBottomDepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
bottomDepth = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'bottomDepth'));
netcdf.close(ncid)
refZBot = -refBottomDepth;
bottomMax = max(refBottomDepth);
gravity = 9.80616;

time = zeros(size(time_fields));
rpeTot = zeros(size(time_fields));
rpeMean = zeros(size(time_fields));
rpePert = zeros(size(time_fields));
keMeanVolume = zeros(size(time_fields));
vertTransportVolume = zeros(size(time_fields));
vertTransportVolumeZ = zeros(nVertLevels,length(time_fields));

fprintf([char(dir) abc ' nt='])
for nt=1:length(time_fields)

  if (multiple_files==true)
    filename = [wd char(dir) abc '/' netcdf_file '_time_' num2str(time_fields(nt)) '.nc'];
    ncid = netcdf.open(filename,'nc_nowrite');
    xtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xtime'));
    hFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'layerThickness'));
    densityFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'potentialDensity'));
    kineticEnergyCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kineticEnergyCell'));
    %vertTransportFull = abs(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'vertTransportVelocityTop')));
    netcdf.close(ncid)
    time_index = 1;  
  else
    time_index = nt;
  end
  
  
% note: time_fields(1) needs to be IC time to compute volume
% correctly.
  if nt==1
  h = squeeze(hFull(:,:,time_index));  
  layerVolume = zeros(nVertLevels,1);
  for iCell=1:nCells
    for k=1:maxLevelCell(iCell)
      layerVolume(k) = layerVolume(k) + h(k,iCell)*areaCell(iCell);
    end
  end 
  layerArea = layerVolume./refLayerThickness;
  end
  
  fprintf([' ' num2str(time_fields(nt))])
  t = xtime(:,time_index)';  % time string
  time(nt)=str2num(t(18:19)) ... % s
         + str2num(t(15:16))   *60 ... %  min
         + str2num(t(12:13))   *60*60 ... %  hour
         +(str2num(t( 9:10))-1)*60*60*24 ... %  day
         +(str2num(t( 6: 7))-1)*60*60*24*30 ... %  month
         + str2num(t( 1: 4))   *60*60*24*30*360;  %  year

  h = squeeze(hFull(:,:,time_index));  
  density = squeeze(densityFull(:,:,time_index));
  
  vol_1D = zeros(1,nVertLevels*nCells);
  density_1D = zeros(1,nVertLevels*nCells);
  zMid_1D = zeros(1,nVertLevels*nCells);
  zMid = zeros(1,nVertLevels);
  j=0;
  for iCell=1:nCells
    kmax = maxLevelCell(iCell);
    zMid(kmax) = -bottomDepth(iCell)+h(kmax,iCell)/2.0;
    for k=kmax-1:-1:1
      zMid(k) = zMid(k+1) + 0.5*(h(k+1,iCell)+h(k,iCell));
    end
    for k=1:maxLevelCell(iCell)
      j = j+1;
      vol_1D(j) = h(k,iCell)*areaCell(iCell);
      density_1D(j) = density(k,iCell);
      zMid_1D(j) = zMid(k);
    end
  end 
  nCells_1D = j;
  totalVolume = sum(vol_1D(1:nCells_1D));
  meanDensity = sum(density_1D(1:nCells_1D).*vol_1D(1:nCells_1D))/totalVolume;
  meanZMid = sum(zMid_1D(1:nCells_1D).*vol_1D(1:nCells_1D))/totalVolume;
  rpeMean(nt) = gravity* meanDensity* meanZMid* totalVolume;

  [sortedDensity,sortedIndex] = sort( density_1D(1:nCells_1D) );
  sortedVolume = vol_1D(sortedIndex);

  sortedRPEtot = zeros(1,nCells_1D);
  sortedRPEpert = zeros(1,nCells_1D);
  sortedZMid = zeros(1,nCells_1D);
  zTop = 0.0;
  k=1; % start at top layer
  for j=1:nCells_1D
    deltaZ = sortedVolume(j)/layerArea(k);
    if zTop-deltaZ<refZBot(k)&k<nVertLevels
      upperVolume = layerArea(k)*(zTop - refZBot(k));
      lowerVolume = sortedVolume(j) - upperVolume;
      deltaZLower = lowerVolume/layerArea(k+1);
      zBot = refZBot(k) - deltaZLower;
      k=k+1;
    else
      zBot = zTop - deltaZ;
    end
    sortedZMid(j) = 0.5*(zTop+zBot);
    %%%%%%%%%%%%%% not sure about bottomMax in next line
    sortedRPEtot(j) = gravity*sortedDensity(j)*(sortedZMid(j)+bottomMax)*sortedVolume(j);
    sortedRPEpert(j) = gravity*(sortedDensity(j)-meanDensity)*(sortedZMid(j)-meanZMid)*sortedVolume(j);
    zTop = zBot;
  end
%figure
%subplot(3,1,1); plot(sortedDensity)
%subplot(3,1,2); plot(sortedZMid+bottomMax)
%subplot(3,1,3); plot(sortedRPEtot)

  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpeTot(nt) = sum(sortedRPEtot)/sum(areaCell);
  rpePert(nt) = sum(sortedRPEpert)/sum(areaCell);

  % note: this only works with no land cells
  keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_index))));
  vertTransportVolume(nt) =0;% mean(mean((vertTransportFull(:,:,time_index))));
  for k=1:nVertLevels
    vertTransportVolumeZ(k,nt) =0;% mean(mean((vertTransportFull(k,:,time_index))));
  end

end
fprintf('\n')

rpeNorm = (rpeTot - rpeTot(1))/rpeTot(1);

DrpeDt = derivative_ord2(time,rpeTot);

n=length(DrpeDt);

meanDrpeDt = mean(DrpeDt(min_n:max_n));
keMeanTime = 0;%mean(keMeanVolume(min_n:max_n));

vertTransportMean = 0;%mean(vertTransportVolume(min_n:max_n));
  for k=1:nVertLevels
    vertTransportMeanZ(k) = 0;%mean((vertTransportVolumeZ(k,:)));
  end
  
