#all: 
#	python driver.py

all: R160km R320km
	
R160km: 
	# R = 160 km
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 160000 --ts 25 --te 30 --dt 2. &> log160km25-30.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 160000 --ts 20 --te 25 --dt 2. &> log160km20-25.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 160000 --ts  5 --te 10 --dt 2. &> log160km5-10.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 160000 --ts  2 --te  7 --dt 2. &> log160km2-7.log & 

R320km:
	# R = 320 km
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 320000 --ts 25 --te 30 --dt 2. &> log320km25-30.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 320000 --ts 20 --te 25 --dt 2. &> log320km20-25.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 320000 --ts  5 --te 10 --dt 2. &> log320km5-10.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 320000 --ts  2 --te  7 --dt 2. &> log320km2-7.log & 

R32km:
	# R = 32 km
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 32000 --ts 25 --te 30 --dt 2. &> log32km25-30.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 32000 --ts 20 --te 25 --dt 2. &> log32km20-25.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 32000 --ts  5 --te 10 --dt 2. &> log32km5-10.log & 
	time python  diffusivity_binning.py -f all_rlzn_particle_data_rads.npz -r 32000 --ts  2 --te  7 --dt 2. &> log32km2-7.log & 

clean:
	trash log*
	trash Case*
