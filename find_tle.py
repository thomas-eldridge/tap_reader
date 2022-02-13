import datetime as dt
import numpy as np
import pdb
from pyorbital.geoloc import ScanGeometry, compute_pixels, get_lonlatalt

def get_tle(time, nimbus):
	""""""
	if nimbus == 'N4':
		tle_file = "nimbus-4.txt"
	elif nimbus == 'N5':
		tle_file = "nimbus-5.txt"
	elif nimbus == 'N6':
		tle_file = "nimbus-6.txt"
	else:
		raise ValueError('Sensor not recognised')
	f = open(tle_file)
	lines = f.readlines()
	lines1 = lines[::2]
	lines2 = lines[1::2]
	dts = []
	for line in lines1:
		epoch_year = int(line[18:20])
		if epoch_year < 50:
			epoch_year += 100
		epoch_day = float(line[20:32])
		true_year = 1900 + epoch_year
		dts.append(dt.datetime(true_year,01,01) + dt.timedelta(days=epoch_day))
	npdts = np.array(dts)
	the_dt = get_dt(time)
	tle_ind = np.where(abs(npdts-the_dt)==np.min(abs(npdts - the_dt)))[0][0]
	tle1 = lines1[tle_ind][:69]
	tle2 = lines2[tle_ind][:69]
	return tle1, tle2
		

def get_dt(time, epoch=dt.datetime(1970,01,01)):
	retval = epoch + dt.timedelta(seconds=time)
	return retval
	
def get_geoloc(time, dpop, nads, roll, pitch, yaw, nimbus, rot=1.25):
	# now works apart from last element (which is nan in demo)
	t_scan_start = (rot/360.)*(180+nads[0])
	t_scan_end = (rot/360.)*(180+nads[-1])
	tle1, tle2 = get_tle(time, nimbus)
	t = get_dt(time)
	scan_points = np.arange(dpop)
	x = np.deg2rad(np.linspace(nads[-1], nads[0], dpop))
	thir = np.vstack((x, np.zeros((len(x),)))).transpose()
	times = np.linspace(t_scan_start, t_scan_end, dpop)
	sgeom = ScanGeometry(thir, times)
	rpy = (roll, pitch, yaw)
	s_times = sgeom.times(t)
	pixels_pos = compute_pixels((tle1,tle2), sgeom, s_times, rpy)
	pos_time = get_lonlatalt(pixels_pos, s_times)
	return pos_time
	
	
