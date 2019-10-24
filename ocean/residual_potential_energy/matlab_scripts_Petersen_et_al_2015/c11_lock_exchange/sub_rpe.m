function [time,rpe,rpeNorm,DrpeDt,meanDrpeDt,keMeanTime,vertTransportMean,vertTransportMeanZ] ...
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
vertTransportFull = abs(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'vertTransportVelocityTop')));
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
vertTransportVolume = zeros(size(time_fields));
vertTransportVolumeZ = zeros(K,length(time_fields));

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

  [density_sorted,sorted_ind] = sort( density_1D(1:nCells_1D) );
  vol_sorted = vol_1D(sorted_ind);

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

  % note: this only works with no land cells
  keMeanVolume(nt) = mean(mean((kineticEnergyCellFull(:,:,time_fields(nt)))));
  vertTransportVolume(nt) = mean(mean((vertTransportFull(:,:, ...
						  time_fields(nt)))));
  for k=1:K
    vertTransportVolumeZ(k,nt) = mean(mean((vertTransportFull(k,:,time_fields(nt)))));
  end
end
rpe(1);

rpeNorm = (rpe - rpe(1))/rpe(1);

DrpeDt = derivative_ord2(time,rpe);

n=length(DrpeDt);
%%%%%%%%%%% make sure this is right
%min_n=1;
%min_n=floor(n/2);
%max_n = floor(n/2);
%max_n = n;

%min_n=90;
%max_n=96;

figure(21)
rho1=1000;
rho2=995;
vel_scale = 1/2*sqrt(9.8*20*(rho1-rho2)/((rho1+rho2)/2)) ;
plot(time/3600,2*sqrt(keMeanVolume),[0 18],vel_scale*[1 1])
xlabel('time, hours')
ylabel('velocity scale, m/s')
legend('from measured KE','theoretical','location','NorthWest')
grid on
set(gcf,'PaperPositionMode','auto','color',[.8 1 .8], ...
	'PaperPosition',[0 0 4.0 3])
fig=[char(dir) abc '_vel_scale_LE_'];
print('-depsc2',['f/' fig '.eps']);
print('-djpeg',['f/' fig '.jpg']);
unix(['epstopdf f/' fig '.eps --outfile=f/' fig '.pdf']);

meanDrpeDt = mean(DrpeDt(min_n:max_n));
keMeanTime = mean(keMeanVolume(min_n:max_n));

vertTransportMean = mean(vertTransportVolume(min_n:max_n));
  for k=1:K
    vertTransportMeanZ(k) = mean((vertTransportVolumeZ(k,:)));
  end

return

subplot(3,1,1); plot(time/3600,rpe); grid on
  subplot(3,1,2); plot(time/3600,rpeNorm); grid on
  subplot(3,1,3); plot(time/3600,DrpeDt,max(time)/3600*[.5 1],meanDrpeDt*[1 1],'--r'); grid on



figure(10)
subplot(2,1,1)
h_3D = reshape(h,[K dims]);
imagesc(squeeze(h_3D(:,2,:)))
colorbar

subplot(2,1,2)
density_3D = reshape(density,[K dims]);
imagesc(squeeze(density_3D(:,2,:)))
colorbar
max(max(max(density)))
min(min(min(density)))

y1=40;
y2=7;
y=1:200;
z0=500;
z1=2000;
fy = z0 + (z1-z0)/2*(1.0+tanh((y-y1)/y2));
figure(2)
subplot(2,1,1)
plot(y,fy)
set(gca,'YDir','reverse')

z = 500:2000;
invfy = y1 + y2*atanh(2*(z-z0)/(z1-z0) - 1) ;
subplot(2,1,2)
plot(invfy,z)
set(gca,'YDir','reverse')
axis([0 200 500 2000])




whos
