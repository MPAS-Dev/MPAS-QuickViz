function [time,rpe,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime] ...
    = sub_rpe( filename, ...
    dims,time_fields,min_n,max_n,title_txt)

% Compute Resting Potential Energy
% Mark Petersen, LANL, Jan 2013

ncid = netcdf.open(filename,'nc_nowrite');

xtime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xtime'));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xCell'));
nCells = size(work,1);
xCell = reshape(work(:,1), dims);

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yCell'));
yCell = reshape(work(:,1), dims);

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'yEdge'));
yEdge = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'xEdge'));
xEdge = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'areaCell'));
areaCell = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'maxLevelCell'));
maxLevelCell = squeeze(work(:,1));

work = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'bottomDepth'));
bottomDepth = squeeze(work(:,1));

hFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'layerThickness'));
densityFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'density'));
kineticEnergyCellFull = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'kineticEnergyCell'));
K = size(hFull,1);

ridgeDepth = 500;
bottomMax = max(bottomDepth);
y1=40e3;
y2=7e3;
yMax = max(yEdge);
yMin = min(yEdge);
xMax = max(xEdge);
gravity = 9.80616;

time = zeros(size(time_fields));
rpe = zeros(size(time_fields));
keMeanVolume = zeros(size(time_fields));

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
  
  % convert volume and density to 1D arrays
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

  % Sort density and volume
  [density_sorted,sorted_ind] = sort( density_1D(1:nCells_1D) );
  vol_sorted = vol_1D(sorted_ind);

  % Integrate RPE from bottom up
  rpe1 = zeros(1,nCells_1D);
  sillWidth = zeros(1,nCells_1D);
  yWidth = zeros(1,nCells_1D);
  z = 0.0;
  for i=1:nCells_1D
    % determine area at z:
    if z>=-ridgeDepth
      sillWidth(i) = yMin;
    else
      sillWidth(i) = max(y1 + y2*atanh(2*(-z-ridgeDepth)/ ...
                  (bottomMax-ridgeDepth) - 1.0) , yMin);
    end
    yWidth(i) = yMax - sillWidth(i);
    area = yWidth(i)*xMax;
    thickness = vol_sorted(i)/area;
    zMid(i) = z - thickness/2;
    z = z - thickness;
    rpe1(i) = gravity*density_sorted(i)*(zMid(i)+bottomMax)*vol_sorted(i);
  end
  % to compare to 2D domain in Ilicak, divide by x-width of domain:
  % divide by surface area to get W/m^2 for dRPE/dt
  rpe(nt) = sum(rpe1)/sum(areaCell);

  keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_fields(nt)))));
end
rpe(1);

rpeNorm = (rpe - rpe(1))/rpe(1);

DrpeDt = derivative_ord2(time,rpe);

n=length(DrpeDt);

meanDrpeDt = mean(DrpeDt(min_n:max_n));
keMeanTime = mean(keMeanVolume(min_n:max_n));


