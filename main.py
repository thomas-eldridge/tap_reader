from Data4to6_new import *
import glob
from netCDF4 import Dataset

def read_TAP_file(filename):
    """Inputs:
        - filename; a string corresponding to the complete path to a Nimbus 4, 5 or 6 TAP file.
    Opens the file in read binary mode. Reads the file and writes it as a NetCDF."""
    if 'Nimbus4' in filename:
        data = Data(filename)
    elif 'Nimbus5' in filename:
        data = Data2(filename)
    elif 'Nimbus6' in filename:
        data = Data2(filename)
    else:
        raise ValueError('Not an N4-6 file')
    return data

def write_NC_file(filename, output_filename=None):
    """Inputs:
        - filename; a string corresponding to the complete path to a Nimbus 4, 5 or 6 TAP file
    Reads the TAP file into a Data object, before writing the output to a NetCDF4 file.
    The NetCDF4 file name will be identical to the TAP file name, but with .TAP replaced by .nc"""
    file_data = read_TAP_file(filename)
    data_fields = Fields(file_data)
    if output_filename == None:
        nc_filename = filename.replace('.TAP','_new.nc')
        nc = Dataset(nc_filename,'w')
    else:
        nc = Dataset(output_filename, 'w')
    nc.mirror_rotation = file_data.od.mirror_rot
    nc.sample_frequency= file_data.od.sample_freq
    nc.orbit_number = file_data.od.orbit_no
    nc.station_code = file_data.od.station_code
    nc.swath_block = file_data.od.swath_block
    nc.swaths_per_record = file_data.od.swaths_per_rec
    nc.locator_number = file_data.od.locator_no
    Y_dim = nc.createDimension('Y', data_fields.data.shape[0])
    X_dim = nc.createDimension('X', data_fields.data.shape[1])
    x_dim = nc.createDimension('x', data_fields.nadangs.shape[1])
    Y_var = create_variable(nc, ['Y'], 'Y', 'i')
    X_var = create_variable(nc, ['X'], 'X', 'i')
    x_var = create_variable(nc, ['x'], 'x', 'i')
    time_var = create_variable(nc, ['Y'], 'time', 'd', 'seconds', 'time since 1970/01/01 00:00:00')
    cell_var = create_variable(nc, ['Y'], 'cell_temp', 'i', 'K', 'detector_cell_temperature')
    elec_var = create_variable(nc, ['Y'], 'electronics_temp', 'i', 'K', 'electronics_temperature')
    refa_var = create_variable(nc, ['Y'], 'ref_temp_A', 'i', 'K', 'housing_temperature')
    refb_var = create_variable(nc, ['Y'], 'ref_temp_B', 'i', 'K', 'housing_temperature')
    refc_var = create_variable(nc, ['Y'], 'ref_temp_C', 'i', 'K', 'housing_temperature')
    refd_var = create_variable(nc, ['Y'], 'ref_temp_D', 'i', 'K', 'housing_temperature')
    roll_var = create_variable(nc, ['Y'], 'roll_error', 'f', 'degrees', 'roll_axis_error')
    pitch_var = create_variable(nc, ['Y'], 'pitch_error', 'f', 'degrees', 'pitch_axis_error')
    yaw_var = create_variable(nc, ['Y'], 'yaw_error', 'f', 'degrees', 'yaw_axis_error')
    height_var = create_variable(nc, ['Y'], 'height', 'i', 'km', 'spacecraft_altitude')
    dpop_var = create_variable(nc, ['Y'], 'data_population', 'i', full_name='scanline pixel number')
    sublat_var = create_variable(nc, ['Y'], 'subsat_lat', 'f', 'degrees_north', 'latitude')
    sublon_var = create_variable(nc, ['Y'], 'subsat_lon', 'f', 'degrees_east', 'longitude')
    flag1_var = create_variable(nc, ['Y'], 'flag_1', 'i', full_name='summary flag: at least one other flag is on')
    flag2_var = create_variable(nc, ['Y'], 'flag_2', 'i', full_name='bad consistency check between sample rate, vehicle time and ground time')
    flag3_var = create_variable(nc, ['Y'], 'flag_3', 'i', full_name='bad vehicle time')
    flag4_var = create_variable(nc, ['Y'], 'flag_4', 'i', full_name='vehicle time inserted by flywheel')
    flag5_var = create_variable(nc, ['Y'], 'flag_5', 'i', full_name='vehicle time carrier is absent')
    flag6_var = create_variable(nc, ['Y'], 'flag_6', 'i', full_name='vehicle time has skipped')
    flag8_var = create_variable(nc, ['Y'], 'flag_8', 'i', full_name='bad sync pulse recognition')
    flag9_var = create_variable(nc, ['Y'], 'flag_9', 'i', full_name='dropout of data signal')
    flag12_var = create_variable(nc, ['Y'], 'flag_12', 'i', full_name='bad swath size')
    anchor_nads_var = create_variable(nc, ['x','Y'], 'anchor_nadang', 'f', 'degrees', 'satellite_viewing_angle_at_anchor_points')
    anchor_lats_var = create_variable(nc, ['x','Y'], 'anchor_lats', 'f', 'degrees_north', 'latitude_of_anchor_points') # remember to add 90
    anchor_lons_var = create_variable(nc, ['x','Y'], 'anchor_lons', 'f', 'degrees_east', 'longitude_of_anchor_points') # remember to change positive direction
    data_var = create_variable(nc, ['X','Y'], 'BBT', 'f', 'K', 'brightness_temperature')
    lats_var = create_variable(nc, ['X','Y'], 'lats_lagrange', 'f', 'degrees_north', 'latitude from interpolation')
    lons_var = create_variable(nc, ['X','Y'], 'lons_lagrange', 'f', 'degrees_east', 'longitude from interpolation')
    lats_var2 = create_variable(nc, ['X','Y'], 'lats_pyorb', 'f', 'degrees_north', 'latitude from pyorbital')
    lons_var2 = create_variable(nc, ['X','Y'], 'lons_pyorb', 'f', 'degrees_east', 'longitude from pyorbital')
    solzen_var = create_variable(nc, ['X','Y'], 'solzen', 'f', 'degrees', 'solar zenith angle')
    satzen_var = create_variable(nc, ['X','Y'], 'satzen', 'f', 'degrees', 'satellite zenith angle')
    solaz_var = create_variable(nc, ['X','Y'], 'solaz', 'f', 'degrees', 'solar azimuth angle')
    sataz_var = create_variable(nc, ['X','Y'], 'sataz', 'f', 'degrees', 'satellite azimuth angle')
    Y_var[:] = np.arange(data_fields.data.shape[0])
    X_var[:] = np.arange(data_fields.data.shape[1])
    time_var[:] = data_fields.truetime
    cell_var[:] = data_fields.cell_temps
    elec_var[:] = data_fields.electro_temps
    refa_var[:] = data_fields.ref_temps_a
    refb_var[:] = data_fields.ref_temps_b
    refc_var[:] = data_fields.ref_temps_c
    refd_var[:] = data_fields.ref_temps_d
    roll_var[:] = data_fields.roll_errors
    pitch_var[:] = data_fields.pitch_errors
    yaw_var[:] = data_fields.yaw_errors
    height_var[:] = data_fields.heights
    dpop_var[:] = data_fields.dpops 
    sublat_var[:] = data_fields.sub_satellite_lats
    sublon_var[:] = data_fields.sub_satellite_lons
    flag1_var[:] = data_fields.flag1
    flag2_var[:] = data_fields.flag2
    flag3_var[:] = data_fields.flag3
    flag4_var[:] = data_fields.flag4
    flag5_var[:] = data_fields.flag5
    flag6_var[:] = data_fields.flag6
    flag8_var[:] = data_fields.flag8
    flag9_var[:] = data_fields.flag9
    flag12_var[:] = data_fields.flag12
    anchor_nads_var[:] = data_fields.nadangs
    anchor_lats_var[:] = data_fields.anchor_lats
    anchor_lons_var[:] = data_fields.anchor_lons
    data_var[:] = data_fields.data
    lats_var[:] = data_fields.lats
    lons_var[:] = data_fields.lons
    lats_var2[:] = data_fields.lats2
    lons_var2[:] = data_fields.lons2
    solzen_var[:] = data_fields.sol_zen
    satzen_var[:] = data_fields.sat_zen
    solaz_var[:] = data_fields.sol_az
    sataz_var[:] = data_fields.sat_az
    nc.close()

def create_variable(dataset, dims, name, dtype, units=None, full_name=None):
    if len(dims)==1:
        var = dataset.createVariable(name, dtype, ((dims[0]),), fill_value=-999)
    elif len(dims)==2:
        var = dataset.createVariable(name, dtype, (dims[1],dims[0],), fill_value=-999)
    else:
        raise ValueError('Not expecting this many dimensions')
    if units!=None:
        var.units = units
    if full_name!=None:
        var.standard_name = full_name
    return var


if __name__ == '__main__':
    write_NC_file('/glusterfs/surft/data/T_Eldridge_data/Nimbus_4_data/window/1970/110/Nimbus4-THIRCH115_1970m0420t003837_o00159_DD15397.TAP')
