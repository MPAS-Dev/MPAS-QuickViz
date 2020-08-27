import os
import subprocess
import yaml
import netCDF4
import numpy as np

########################################################################
########################################################################

def write_coordinate_file(grid_file):

  # Read in mesh
  grid_file = grid_file
  grid_nc = netCDF4.Dataset(grid_file,'r')
  lon_grid = np.degrees(grid_nc.variables['lonCell'][:])
  lat_grid = np.degrees(grid_nc.variables['latCell'][:])
  nCells = len(grid_nc.dimensions['nCells'])
  print(nCells)

  # Write coordinate file for OTPS2
  f = open('lat_lon','w')
  for i in range(nCells):
    f.write(str(lat_grid[i])+'  '+str(lon_grid[i])+'\n')
  f.close()

########################################################################
########################################################################

def setup_otps2(constituents,tpxo8_data_path):

  for con in constituents:
    print('setup '+con)
 
    # Lines for the setup_con files
    lines = [{'inp':'inputs/Model_atlas_'+con ,      'comment':'! 1. tidal model control file'},
             {'inp':'lat_lon' ,                      'comment':'! 2. latitude/longitude/<time> file'},
             {'inp':'z' ,                            'comment':'! 3. z/U/V/u/v'},
             {'inp': con ,                           'comment':'! 4. tidal constituents to include'},
             {'inp':'AP' ,                           'comment':'! 5. AP/RI'},
             {'inp':'oce' ,                          'comment':'! 6. oce/geo'},
             {'inp':'1' ,                            'comment':'! 7. 1/0 correct for minor constituents'},
             {'inp':'outputs/'+con+'.out' ,          'comment':'! 8. output file (ASCII)'}]
    
    # Create directory for setup_con and Model_atlas_con files
    if not os.path.exists('inputs'):
      os.mkdir('inputs') 
   
    # Write the setup_con file
    f = open('inputs/'+con+'_setup','w')
    for line in lines:
      spaces = 28 - len(line['inp'])
      f.write(line['inp'] + spaces*' ' + line['comment'] + '\n')
    f.close()

    # Write the Model_atlas_con file
    f = open('inputs/Model_atlas_'+con,'w')
    f.write('TPXO8/hf.'+con+'_tpxo8_atlas_30c_v1.out\n')
    f.write('TPXO8/uv.'+con+'_tpxo8_atlas_30c_v1.out\n')
    f.write('TPXO8/grid_tpxo8_atlas_30_v1')
    f.close()

    # Link the TPXO8 data directory
    subprocess.call('ln -sf ' + tpxo8_data_path + ' TPXO8', shell=True)

    # Create directory for the con.out files
    if not os.path.exists('outputs'):
      os.mkdir('outputs') 

########################################################################
########################################################################

def run_otps2(exe_path,constituents):

  # Make the executable if necessary 
  if not os.path.isfile(exe_path+'/extract_HC'):
    pwd = os.getcwd()
    os.chdir(exe_path)
    subprocess.call('make extract_HC',shell=True)
    os.chdir(pwd)

  # Run the executable 
  for con in constituents:
    print('run '+con)
    subprocess.call(exe_path+'/extract_HC < inputs/'+con+'_setup',shell=True)

########################################################################
########################################################################

def read_otps2_output(constituents):

  bou_AP = {}
  for con in constituents:
    bou_AP[con] = {'amp':[], 'phase':[]}

    f = open('outputs/'+con+'.out','r')
    lines = f.read().splitlines()
    for line in lines[3:]:
      line_sp = line.split()
      if line_sp[2] != '*************':
        val = float(line_sp[2])
        bou_AP[con]['amp'].append(val)
      else:
        bou_AP[con]['amp'].append('-9999')

      if line_sp[3] != 'Site':
        val = float(line_sp[3])
        if val < 0:
          val = val + 360.0
        bou_AP[con]['phase'].append(val)
        
      else:
        bou_AP[con]['phase'].append(-9999)

    #pprint.pprint(bou_AP)

  return bou_AP

########################################################################
########################################################################

def append_tpxo8_data(output_file,constituents,mesh_AP):

  data_nc = netCDF4.Dataset(output_file,'a', format='NETCDF3_64BIT_OFFSET')
  for con in constituents:
    amp_var = data_nc.createVariable(con.upper()+'AmplitudeTPXO8',np.float64,('nCells'))
    amp_var[:] = mesh_AP[con]['amp'][:]
    amp_var.units = 'm'
    amp_var.long_name = 'Amplitude of '+con.upper()+ ' tidal consitiuent at each cell center from TPXO8 model'
    phase_var = data_nc.createVariable(con.upper()+'PhaseTPXO8',np.float64,('nCells'))
    phase_var[:] = mesh_AP[con]['phase'][:]
    phase_var.units = 'deg'
    phase_var.long_name = 'Phase of '+con.upper()+ ' tidal consitiuent at each cell center from TPXO8 model'


########################################################################
########################################################################

if __name__ == '__main__':

  pwd = os.getcwd()
  inputfile = pwd+'/inject_TPXO8.config'
  f = open(inputfile)
  cfg = yaml.load(f,yaml.Loader)

  write_coordinate_file(cfg['grid_file'])  
  setup_otps2(cfg['constituents'],cfg['tpxo8_data_path'])
  run_otps2(cfg['otps2_exe_path'],cfg['constituents'])
  mesh_AP = read_otps2_output(cfg['constituents'])
  append_tpxo8_data(cfg['output_file'],cfg['constituents'],mesh_AP)

