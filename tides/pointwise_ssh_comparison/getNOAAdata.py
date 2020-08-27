import csv
import pprint
import urllib2
import calendar
from bs4 import BeautifulSoup
import datetime
import json
import subprocess

station_html_file = './tidesandcurrents.noaa.gov_stations.html'
product_ID = 'water_level'
#product_ID = 'currents'
organization = 'LANL'
year = 2012
months = [10,11]
states = ['NY','NJ','DE','MD','VA']

###################################################################################################
###################################################################################################

def read_station_html_file(station_html_file,states):
  station_dict = {}
  f = open(station_html_file)
  soup = BeautifulSoup(f,features="lxml")
  states = soup.find_all(attrs={'class':'span12 areaheader'})    
  for state in states:
    stations = state.find_all(attrs={'class':'span4'})
    for station in stations:
      sta = station.find('a').text.split()
      dates = station.find('span').text.split("-") 
      
      st = sta[-1]
      if st not in station_dict:
        station_dict[st] = []
      
      n = len(sta)
      ID = sta[0]
      name = ' '.join(sta[1:n])
      start = dates[0].strip()
      end = dates[1].strip()
      
      station_dict[st].append({'id':ID,'name':name,'start':start,'end':end})    
      
  pprint.pprint(station_dict)    
  return station_dict

###################################################################################################
###################################################################################################

if __name__ == '__main__':

 
  # API parameters for data retrieval
  base_url = 'https://tidesandcurrents.noaa.gov/api/datagetter?'
  product = 'product=' + product_ID 
  datum = 'datum=' + 'MSL'
  units = 'units=' + 'metric'
  time_zone = 'time_zone=' + 'gmt'
  application = 'application=' + organization 
  frmt = 'format=' + 'json'
 
  # Read in NOAA station file
  station_dict = read_station_html_file(station_html_file,states) 
  
  # Determine requested start and end dates
  sim_start = datetime.datetime(year,months[0],1)
  sim_end = datetime.datetime(year,months[-1],calendar.monthrange(year,months[-1])[1])
  
  station_file = open('stations.txt','w+')
  
  for st in states:
  
    stations = station_dict[st]
    
    for station in stations:        
       
      if not station['start']:
        sta_start = datetime.datetime(1950,1,1)
      else:
        sta_start = datetime.datetime.strptime(station['start'],'%b %d, %Y')
        
      if station['end'] == 'present':
        sta_end = datetime.datetime.now()
      else:  
        sta_end =  datetime.datetime.strptime(station['end'],'%b %d, %Y')
      
      # Decide if station has data in the requested time-frame
      if sta_start < sim_start and sta_end > sim_end:
    
        ID = station['id']
        error = False
        print station

        observation_filename = ID+'_'+str(year)+'.txt'          
        observation_file = open(observation_filename,'w+')    
        observation_file.write('YYYY MM DD hh mm SSH\n')

        for month in months:
          # Create additional API parameters 
          sta_ID = 'station='+station['id']     
          m = str(month).zfill(2)   
          begin_date = 'begin_date='+str(year)+m+'01'
          d = str(calendar.monthrange(year,month)[1]).zfill(2)
          end_date = 'end_date='+str(year)+m+d
          
          # Build full URL 
          info = '&'.join([begin_date,end_date,sta_ID,product,datum,units,time_zone,application,frmt])
          url = base_url + info

          # Download and read in data
          data = urllib2.urlopen(url).read()
          data_dict = json.loads(data)
          
          # Check for errors
          if 'error' in data_dict:
            print data_dict['error']['message']
            error = True
          else: 
            # Get station metadata
            lat = data_dict['metadata']['lat']
            lon = data_dict['metadata']['lon']
            name = data_dict['metadata']['name']                
            
            # Get and write station data
            for line in data_dict['data']:
              date = line['t']
              obs = line['v']
              if not obs:
                obs = "99.0"
              dt = datetime.datetime.strptime(date,'%Y-%m-%d %H:%M') 
              t  = datetime.datetime.strftime(dt,'%Y %m %d %H %M')
              observation_file.write(t+" "+obs+"\n")           

        observation_file.close()
        
        # Write to station list file
        if not error:
          station_file.write(' '.join([lon,lat,ID,name,'\n']))
        else:
          subprocess.call(['rm',observation_filename])        
