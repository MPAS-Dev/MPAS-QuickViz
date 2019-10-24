function [time,rpeTot,rpePert,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime] ...
    = sub_rpe( ...
    wd,dir,abc,hist_files, kmt_file, grid_file, depth_file, ...
    dims,time_fields,min_n,max_n,title_txt)


linear_eos = false;
gravity = 9.80616;

% Compute Resting Potential Energy
% Mark Petersen, LANL, Jan 2013

px = [.53];
py=[.52]; % Midpoint position of plots
pw = [.87];  ph=[.75]; % width and height of plots

filename = [wd char(dir) abc '/' char(hist_files(1))];
ncid = netcdf.open(filename,'nc_nowrite');
[dimname, nx] = netcdf.inqDim(ncid,0);
[dimname, ny] = netcdf.inqDim(ncid,1);
[dimname, nz] = netcdf.inqDim(ncid,2);
netcdf.close(ncid)

depth_data = load([wd char(dir) abc '/' depth_file]);
dz = depth_data(:,1)/100; % refLayerThickness in m
bottomDepth = depth_data(:,3); % z coord at bottom of layer in m
refZBot = -depth_data(:,3); % z coord at bottom of layer in m
maxBottom = sum(dz);

fid = fopen([wd char(dir) abc '/' kmt_file]);
  KMT = fread(fid,[nx, ny], 'int','ieee-be');
fclose(fid);
maxLevelCell = KMT; % convert from POP variable to mpas variable

fid = fopen([wd char(dir) abc '/' grid_file]);
   ULAT = fread(fid,[nx, ny], 'real*8','ieee-be');
   ULONG = fread(fid,[nx, ny], 'real*8','ieee-be');
   HTN = fread(fid,[nx, ny], 'real*8','ieee-be')/100.0; % in m
   HTE = fread(fid,[nx, ny], 'real*8','ieee-be')/100.0; % in m
   HUS = fread(fid,[nx, ny], 'real*8','ieee-be')/100.0; % in m
   HUW = fread(fid,[nx, ny], 'real*8','ieee-be')/100.0; % in m
   ANGLE = fread(fid,[nx, ny], 'real*8','ieee-be');
fclose(fid);

DXU = zeros(nx,ny);
DYU = zeros(nx,ny);
DXT = zeros(nx,ny);
DYT = zeros(nx,ny);
UAREA = zeros(nx,ny);
TAREA = zeros(nx,ny);
for i=1:nx-1
  for j=1:ny
    DXU(i,j) =  0.5*(HTN(i,j)+HTN(i+1,j));
  end
end
DXU(nx,:) = DXU(nx-1,:);
for i=1:nx
  for j=1:ny-1
    DYU(i,j) =  0.5*(HTE(i,j)+HTE(i,j+1));
  end
end
DYU(:,ny) = DYU(:,ny-1);
for i=1:nx
  for j=2:ny
    DXT(i,j) =  0.5*(HTN(i,j)+HTN(i,j-1));
  end
end
DXT(:,1) = DXT(:,2);
for i=2:nx
  for j=1:ny
    DYT(i,j) =  0.5*(HTE(i,j)+HTE(i-1,j));
  end
end
DYT(1,:) = DYT(2,:);
UAREA = DXU.*DYU;
TAREA = DXT.*DYT;
areaCell = TAREA;  % convert from POP variable to mpas variable,
                   % cell area in m^2

keMeanVolume = zeros(1,length(hist_files));

  ncid = netcdf.open([wd char(dir) abc '/' char(hist_files(1))],'nc_nowrite');
    PD = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'PD'))*1000; % in kg/m^3
    %PD2000 = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'PD2000'));
    SSH = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SSH'))*1000; % in kg/m^3
  netcdf.close(ncid)

% note: time_fields(1) needs to be IC time to compute volume correctly.
  layerVolume = zeros(nz,1);
  for i=1:nx
    for j=1:ny
      for k=1:KMT(i,j)
	% note: for z-star, need to alter this by the ssh.
        layerVolume(k) = layerVolume(k) + dz(k)*TAREA(i,j);
       end
    end 
  end 
layerArea = layerVolume./dz;

fprintf(['begin, ' char(dir) abc ' nt='])
for nt=1:length(time_fields)
  fprintf([' ' num2str(nt)])

  if (linear_eos)
        config_eos_linear_alpha = 0.2;
        config_eos_linear_beta = 0.8;
        config_eos_linear_Tref = 5.0;
        config_eos_linear_Sref = 35.0;
        config_eos_linear_densityref = 1000.0;
     density = config_eos_linear_densityref - ...
	       config_eos_linear_alpha*(TEMP-config_eos_linear_Tref);
  else
    density = squeeze(PD(:,:,:,nt));
    ind = find(density<0);
    density(ind)=NaN;
  end

  id=0;
  h = zeros(nx,ny,nz);
  vol_1D = zeros(1,nx*ny*nz);
  density_1D = zeros(1,nx*ny*nz);
  zMid_1D = zeros(1,nx*ny*nz);
  zMid = zeros(1,nz);
  for i=1:nx
    for j=1:ny
      if KMT(i,j)>0
      for k=1:KMT(i,j)
	% note: for z-star, need to alter this by the ssh.
	h(i,j,k) = dz(k);
      end
      kmax = KMT(i,j);
      zMid(kmax) = -bottomDepth(kmax)+h(i,j,k)/2.0;
      for k=kmax-1:-1:1
        zMid(k) = zMid(k+1) + 0.5*(h(i,j,k+1)+h(i,j,k));
      end
      for k=1:KMT(i,j)
        id = id+1;
        vol_1D(id) = h(i,j,k)*TAREA(i,j); % in m^3
        density_1D(id) = density(i,j,k); % in kg/m^3
        zMid_1D(id) = zMid(k);
      end
      end
    end 
  end 
  nCells_1D = id;
  totalVolume = sum(vol_1D(1:nCells_1D));
  meanDensity = sum(density_1D(1:nCells_1D).*vol_1D(1:nCells_1D))/totalVolume;
  meanZMid = sum(zMid_1D(1:nCells_1D).*vol_1D(1:nCells_1D))/totalVolume;
  rpeMean(nt) = gravity* meanDensity* meanZMid* totalVolume;

  [sortedDensity,sorted_ind] = sort( density_1D(1:nCells_1D) );
  sortedVolume = vol_1D(sorted_ind);

  sortedRPEtot = zeros(1,nCells_1D);
  sortedRPEpert = zeros(1,nCells_1D);
  sortedZMid = zeros(1,nCells_1D);
  zTop = 0.0;
  k=1; % start at top layer
  for j=1:nCells_1D
    deltaZ = sortedVolume(j)/layerArea(k);
    if zTop-deltaZ<refZBot(k)&k<nz
      upperVolume = layerArea(k)*(zTop - refZBot(k));
      lowerVolume = sortedVolume(j) - upperVolume;
      deltaZLower = lowerVolume/layerArea(k+1);
      zBot = refZBot(k) - deltaZLower;
      k=k+1;
    else
      zBot = zTop - deltaZ;
    end
    sortedZMid(j) = 0.5*(zTop+zBot);
    %%%%%%%%%%%%%% not sure about bottomDepth(nz) in next line
    sortedRPEtot(j)  = gravity*sortedDensity(j)*(sortedZMid(j)+bottomDepth(nz))*sortedVolume(j);
    sortedRPEpert(j) = gravity*(sortedDensity(j)-meanDensity)*(sortedZMid(j)-meanZMid)*sortedVolume(j);
    zTop = zBot;
  end
  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpeTot(nt) = sum(sortedRPEtot)/layerArea(1);
  rpePert(nt) = sum(sortedRPEpert)/layerArea(1);

  
if (1==2) % this is for flat bottom only
  zMid = zeros(1,nCells_1D);
  rpeSorted = zeros(1,nCells_1D);
  % z is positive in this calculation
  z = maxBottom;
  % Here sumCellArea is a constant, so
  % this only works for flat-bottom cases.  For topography, see m55_global_spindown
  sumCellArea = sum(sum(min(KMT,1).*areaCell));
  for i=1:nCells_1D
    thickness = vol_sorted(i)/sumCellArea;
    zMid(i) = z - thickness/2;
    z = z - thickness;
    rpeSorted(i) = gravity*sortedDensity(i)*zMid(i)*vol_sorted(i);
  end
  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpeTot(nt) = sum(rpeSorted)/sumCellArea;
end % this is for flat bottom only
  
%  volSum=0;
%  keSum=0;
%  for i=1:nx
%    for j=2:ny
%      for k=1:KMT(i,j)
%	vol = dz(k)*UAREA(i,j);
%        volSum = volSum + vol;
%	keSum = keSum + (UVEL(i,j,k)^2 + VVEL(i,j,k)^2)*vol;
%      end
%    end 
%  end 
%  keMeanVolume(nt) = keSum/2.0/volSum;
  
end

% this is a hack: time in seconds, originally 5 days.
time = 5*[1:length(time_fields)] *60*60*24
rpeTot
rpeNorm = (rpeTot - rpeTot(1))/rpeTot(1)
DrpeDt = derivative_ord2(time,rpeTot)
meanDrpeDt = mean(DrpeDt(min_n:max_n))

%keMeanTime = mean(keMeanVolume(min_n:max_n));
keMeanTime = 0.0;

return


%   % note: this only works with no land cells
%   keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_fields(nt)))));
% %  vertTransportVolume(nt) = mean(mean((vertTransportFull(:,:, ...
% %						  time_fields(nt)))));
% %  for k=1:K
% %    vertTransportVolumeZ(k,nt) = mean(mean((vertTransportFull(k,:,time_fields(nt)))));
% %  end

% areaCell = squeeze(work(:,1));

% work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'maxLevelCell'));
% maxLevelCell = squeeze(work(:,1));

% work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'refBottomDepth'));
% refBottomDepth = squeeze(work(:,1));

% hFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'layerThickness'));
% densityFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'potentialDensity'));
% kineticEnergyCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kineticEnergyCell'));
% %vertTransportFull = abs(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'vertTransportVelocityTop')));
% K = size(hFull,1);
% nCells = size(areaCell,1);
% maxBottom = max(refBottomDepth);
% gravity = 9.80616;

% time = zeros(size(time_fields));
% rpe = zeros(size(time_fields));
% keMeanVolume = zeros(size(time_fields));
% %vertTransportVolume = zeros(size(time_fields));
% %vertTransportVolumeZ = zeros(K,length(time_fields));

% for nt=1:length(time_fields)
%   t = xtime(:,time_fields(nt))';  % time string
%   time(nt)=str2num(t(18:19)) ... % s
%          + str2num(t(15:16))   *60 ... %  min
%          + str2num(t(12:13))   *60*60 ... %  hour
%          +(str2num(t( 9:10))-1)*60*60*24 ... %  day
%          +(str2num(t( 6: 7))-1)*60*60*24*30 ... %  month
%          + str2num(t( 1: 4))   *60*60*24*30*360;  %  year

%   h = squeeze(hFull(:,:,time_fields(nt)));  
%   density = squeeze(densityFull(:,:,time_fields(nt)));
  
%   vol_1D = zeros(1,K*nCells);
%   density_1D = zeros(1,K*nCells);
%   i=0;
%   for iCell=1:nCells
%     for k=1:maxLevelCell(iCell)
%       i = i+1;
%       vol_1D(i) = h(k,iCell)*areaCell(iCell);
%       density_1D(i) = density(k,iCell);
%     end
%   end 
%   nCells_1D = i;

% %  fprintf('sorting...')
%   [sortedDensity,sorted_ind] = sort( density_1D(1:nCells_1D) );
% %  fprintf('done \n')
%   vol_sorted = vol_1D(sorted_ind);

%   zMid = zeros(1,nCells_1D);
%   rpeSorted = zeros(1,nCells_1D);
%   z = maxBottom;
%   % this only works for flat-bottom cases:
%   sumCellArea = sum(areaCell);
%   for i=1:nCells_1D
%     thickness = vol_sorted(i)/sumCellArea;
%     zMid(i) = z - thickness/2;
%     z = z - thickness;
%     rpeSorted(i) = gravity*sortedDensity(i)*zMid(i)*vol_sorted(i);
%   end
%   % to compare to 2D domain in Ilicak, divide by x-width of domain:
%   % divide by surface area to get W/m^2 for dRPE/dt
%   rpe(nt) = sum(rpeSorted)/sum(areaCell);

%   % note: this only works with no land cells
%   keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_fields(nt)))));
% %  vertTransportVolume(nt) = mean(mean((vertTransportFull(:,:, ...
% %						  time_fields(nt)))));
% %  for k=1:K
% %    vertTransportVolumeZ(k,nt) = mean(mean((vertTransportFull(k,:,time_fields(nt)))));
% %  end
% end
% rpe(1);

% rpeNorm = (rpeTot - rpeTot(1))/rpeTot(1);
% DrpeDt = derivative_ord2(time,rpeTot);
% meanDrpeDt = mean(DrpeDt(min_n:max_n));

% n=length(DrpeDt);
% %%%%%%%%%%% make sure this is right
% %min_n=1;
% %min_n=floor(n/2);
% %max_n = floor(n/2);
% %max_n = n;

% %min_n=90;
% %max_n=96;

% meanDrpeDt = mean(DrpeDt(min_n:max_n));
% keMeanTime = mean(keMeanVolume(min_n:max_n));

% %vertTransportMean = mean(vertTransportVolume(min_n:max_n));
% %  for k=1:K
% %    vertTransportMeanZ(k) = mean((vertTransportVolumeZ(k,:)));
% %  end

% vertTransportMean=0;
% vertTransportMeanZ=0;
