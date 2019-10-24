function [time,rpeTot,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime,keMeanVolume] ...
    = sub_rpe( ...
    wd,dir,abc,hist_files, kmt_file, grid_file, depth_file, ...
    time_fields,min_n,max_n,title_txt)

linear_eos = true;
gravity = 9.80616;

% Compute Resting Potential Energy
% Mark Petersen, LANL, Jan 2013

px = [.53];
py=[.52]; % Midpoint position of plots
pw = [.87];  ph=[.75]; % width and height of plots

filename = [wd char(dir) abc '/h.001.nc']
ncid = netcdf.open(filename,'nc_nowrite');
[dimname, nx] = netcdf.inqDim(ncid,0);
[dimname, ny] = netcdf.inqDim(ncid,1);
[dimname, nz] = netcdf.inqDim(ncid,2);
netcdf.close(ncid)

depth_data = load(['data/' depth_file]);
dz = depth_data(:,1)/100; % refLayerThickness in m
maxBottom = sum(dz);

fid = fopen(['data/' kmt_file]);
  KMT = fread(fid,[nx, ny], 'int','ieee-be');
fclose(fid);
maxLevelCell = KMT; % convert from POP variable to mpas variable

fid = fopen(['data/' grid_file]);
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
fprintf('computing volume...')
vol = zeros(nx,ny,nz);
  for i=1:nx
    for j=1:ny
      for k=1:KMT(i,j)
	vol(i,j,k) = UAREA(i,j)*dz(k);
      end
    end 
  end 
  volSum = sum(sum(sum(vol)));
fprintf('done \n')
		   
keMeanVolume = zeros(1,length(time_fields));

for nt=1:length(time_fields)
  fprintf('begin, nt=\n')
nt
  ncid = netcdf.open([wd char(dir) abc '/h.' num2str_fixed0(nt,'%g',3) '.nc'],'nc_nowrite');
    TEMP = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'TEMP'));
    UVEL = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'UVEL'))/100.0; %m/s
    VVEL = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'VVEL'))/100.0; %m/s
    SHGT = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'SHGT'));
  netcdf.close(ncid)
  fprintf('after read TEMP,UVEL,VVEL \n')

  if (linear_eos)
        config_eos_linear_alpha = 0.2;
        config_eos_linear_beta = 0.8;
        config_eos_linear_Tref = 5.0;
        config_eos_linear_Sref = 35.0;
        config_eos_linear_densityref = 1000.0;
     density = config_eos_linear_densityref - config_eos_linear_alpha*(TEMP-config_eos_linear_Tref);
  end

  id=0;
  h = zeros(nx,ny,nz);
  vol_1D = zeros(1,nx*ny*nz);
  density_1D = zeros(1,nx*ny*nz);
  for i=1:nx
    for j=1:ny
      for k=1:KMT(i,j)
        id = id+1;
	% note: for z-star, need to alter this by the ssh.
	h(i,j,k) = dz(k);
        vol_1D(id) = h(i,j,k)*TAREA(i,j); % in m^3
        density_1D(id) = density(i,j,k); % in kg/m^3
      end
    end 
  end 
  nCells_1D = id;

  [density_sorted,sorted_ind] = sort( density_1D(1:nCells_1D) );
  vol_sorted = vol_1D(sorted_ind);
  fprintf('after sort \n')

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
    rpeSorted(i) = gravity*density_sorted(i)*zMid(i)*vol_sorted(i);
  end
  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpeTot(nt) = sum(rpeSorted)/sumCellArea;

 fprintf('after RPE calc \n')

  ke = (UVEL.^2 + VVEL.^2).*vol;
  keSum=sum(sum(sum(ke)));
  keMeanVolume(nt) = keSum/2.0/volSum;

  fprintf('after KE calc \n')

end

% time in seconds
time = time_fields*60*60*24
rpeTot
rpeNorm = (rpeTot - rpeTot(1))/rpeTot(1)
DrpeDt = derivative_ord2(time,rpeTot)
meanDrpeDt = mean(DrpeDt(min_n:max_n))

keMeanTime = mean(keMeanVolume(min_n:max_n));

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
%   [density_sorted,sorted_ind] = sort( density_1D(1:nCells_1D) );
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
%     rpeSorted(i) = gravity*density_sorted(i)*zMid(i)*vol_sorted(i);
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
