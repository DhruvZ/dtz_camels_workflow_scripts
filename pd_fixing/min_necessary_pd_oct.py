import numpy as np
import yt
from hyperion.model import Model
from astropy import units as u
from astropy.modeling.models import BlackBody
import h5py

galaxy_num = 569
fname = f'/orange/narayanan/d.zimmerman/camels_test/filtered/SIMBA_1P_2_0/snap33/galaxy_{galaxy_num}.hdf5'# name of filtered file
nref = 16
dustdir = '/home/desika.narayanan/hyperion-dust-0.1.0/dust_files/' #location of your dust files
dustfile = 'd03_3.1_6.0_A.hdf5'
bbox_lim = 1.e5
x_cent = 19496.465
y_cent = 9300.578
z_cent = 6220.7397
TCMB = 2.7300000000000004
model = f'/orange/narayanan/d.zimmerman/camels_test/pd_scripts/SIMBA_1P_2_0/snap33/snap33_galaxy{galaxy_num}.py'


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()

    return idx

def get_J_CMB():
    #returns the mean intensity for the CMB integrated over min_lam to
    #max_lam (i.e. returns erg/s/cm**2; the same thing as doing 4*sigma
    #T^4)
    min_lam = (1.*u.angstrom).to(u.micron)
    max_lam = (1*u.cm).to(u.micron)

    wavelengths = np.linspace(min_lam,max_lam,1.e5)

    bb_lam = BlackBody(TCMB * u.K, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
    flux = bb_lam(wavelengths)

    J = np.trapz(flux,wavelengths).to(u.erg/u.s/u.cm**2/u.sr)
    solid_angle = 4.*np.pi*u.sr
    J = J*solid_angle
    return J

def energy_density_absorbed_by_CMB():
    extinction_file = dustdir+dustfile

    mw_df = h5py.File(extinction_file,'r')
    mw_o = mw_df['optical_properties']
    mw_df_nu = mw_o['nu']*u.Hz
    mw_df_chi = mw_o['chi']*u.cm**2/u.g
    bb_nu = BlackBody(TCMB * u.K)
    b_nu = bb_nu(mw_df_nu)
    #energy_density_absorbed = 4pi int b_nu * kappa_nu d_nu since b_nu
    #has units erg/s/cm^2/Hz/str and kappa_nu has units cm^2/g.  this results in units erg/s/g
    steradians = 4*np.pi*u.sr
    energy_density_absorbed = (steradians*np.trapz( (b_nu*mw_df_chi),mw_df_nu)).to(u.erg/u.s/u.g)
    return energy_density_absorbed

def field_add(fname, bounding_box=None, ds=None,add_smoothed_quantities=True):

    def _gassmoothinglength(field,data):
        return data[('PartType0', 'SmoothingLength')].in_units('pc')

    def _starmetals_00(field, data):
        el_dict = {'He': '01',
                   'C': '02',
                   'N': '03',
                   'O': '04',
                   'Ne': '05',
                   'Mg': '06',
                   'Si': '07',
                   'S': '08',
                   'Ca': '09',
                   'Fe': '10'}
        el_str = field.name[1]
        if '_' in el_str:
            el_name = field.name[1][field.name[1].find('_')+1:]
            el_num = el_dict[el_name]
        else:
            el_num = '00'
        return data[('PartType4', 'Metallicity_'+el_num)]
    def _starmetals(field, data):
        return data[('PartType4', 'Metallicity')]

    def _starcoordinates(field, data):
        return data[('PartType4', 'Coordinates')]

    def _starformationtime(field, data):
        return data[('PartType4', 'StellarFormationTime')]

    def _starmasses(field, data):
        return data[("PartType4", "Masses")]

    def _gasdensity(field, data):
        return data[('PartType0', 'Density')]
    def _gasmetals_00(field, data):
        el_dict = {'He': '01',
                   'C': '02',
                   'N': '03',
                   'O': '04',
                   'Ne': '05',
                   'Mg': '06',
                   'Si': '07',
                   'S': '08',
                   'Ca': '09',
                   'Fe': '10'}
        el_str = field.name[1]
        if '_' in el_str:
            el_name = field.name[1][field.name[1].find('_')+1:]
            el_num = el_dict[el_name]
        else:
            el_num = '00'

        return data[('PartType0', 'Metallicity_'+el_num)]

    def _gasmetals(field, data):
        return data[('PartType0', 'Metallicity')]

    def _gascoordinates(field, data):
        return data[('PartType0', 'Coordinates')]

    def _gasmasses(field, data):
        return data[('PartType0', 'Masses')]

    def _gasfh2(field, data):
        try: return data[('PartType0', 'FractionH2')]
        except: return data[('PartType0', 'metallicity')]*0.

    def _gassfr(field, data):
        return data[('PartType0', 'StarFormationRate')]

    def _gassmootheddensity(field, data):
        return data.ds.parameters['octree'][('PartType0', 'density')]

    def _gassmoothedmetals(field, data):        
        try:
            el_str = field.name[1]
            if '_' in el_str:
                el_name = field.name[1][field.name[1].find('_')+1:]+"_"
            else:
                el_name = ""
            return data.ds.parameters['octree'][('PartType0', el_name+'metallicity')]
        except:
            return data.ds.parameters['octree'][('PartType0', 'metallicity')]

    def _gassmoothedmasses(field, data):
        if float(yt.__version__[0:3]) >= 4:
            return data.ds.parameters['octree'][('PartType0', 'Masses')]
        else:
            return data[('deposit', 'PartType0_mass')]

    def _metaldens_00(field, data):
        return (data["PartType0", "Density"]*data["PartType0", "Metallicity_00"])

    def _metaldens(field, data):
        return (data["PartType0", "Density"]*data["PartType0", "Metallicity"])

    def _metalmass_00(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "Metallicity_00"].value))

    def _metalmass(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "Metallicity"].value))

    def _metalsmoothedmasses(field, data):
        return (data.ds.parameters['octree'][('PartType0', 'Masses')]* data.ds.parameters['octree'][('PartType0','metallicity')])

    def _dustmass_manual(field, data):
        return (data.ds.arr(data[("PartType0", "Dust_Masses")].value, 'code_mass'))

    def _dustsmoothedmasses(field, data):
        dsm = ds.arr(data.ds.parameters['octree'][('PartType0','Dust_Masses')],'code_mass')
        return dsm

    def _return_dust_mass(field,data):
        return data['dust','mass']
    def _stellarages(field, data):
        ad = data.ds.all_data()
        if data.ds.cosmological_simulation == False:
            simtime = data.ds.current_time.in_units('Gyr')
            simtime = simtime.value

            age = simtime-data.ds.arr(ad[('PartType4', 'StellarFormationTime')],'Gyr').value
            # make the minimum age 1 million years
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print(
                '[gadget2pd: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ', simtime)
            print('--------------\n')
        else:
            yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=data.ds.hubble_constant,
                                                        omega_matter=data.ds.omega_matter,
                                                        omega_lambda=data.ds.omega_lambda)
            simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units('Gyr').value # Current age of the universe
            scalefactor = data[('PartType4', 'StellarFormationTime')].value
            formation_z = (1./scalefactor)-1.
            formation_time = yt_cosmo.t_from_z(formation_z).in_units('Gyr').value
            age = simtime - formation_time
            # Minimum age is set to 1 Myr (FSPS doesn't work properly for ages below 1 Myr)
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print('[gadget2pd: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Gyr): ', simtime)

        age = data.ds.arr(age, 'Gyr')
        return age

    def _starsmoothedmasses(field, data):
        return data[('deposit', 'PartType4_mass')]


    def _size_with_units(field,data):
        return data.ds.parameters['size']


    # load the ds (but only if this is our first passthrough and we pass in fname)
    if fname != None:
        ds = yt.load(fname, bounding_box=bounding_box)
        
        #ds.sph_smoothing_style = "gather"
        ds.index
        ad = ds.all_data()
        #ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)

    #if we're in the 4.x branch of yt, load up the octree for smoothing
    left = np.array([pos[0] for pos in bounding_box])
    right = np.array([pos[1] for pos in bounding_box])
    #octree = ds.octree(left, right, over_refine_factor=cfg.par.oref, n_ref=cfg.par.n_ref, force_build=True)
    octree = ds.octree(left,right,n_ref=nref)
    ds.parameters['octree'] = octree

    print ('BOUNDING BOX:', bounding_box, 'LEFT: ', left, 'RIGHT: ', right)

    # for the metal fields have a few options since gadget can have different nomenclatures
    ad = ds.all_data()

    if ('PartType4', 'Metallicity_00') in ds.derived_field_list:
        try:
            ds.add_field(('star','metals'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_He'), function=_starmetals_00, sampling_type='particle', units="code_metallicity")
            ds.add_field(('star','metals_C'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_N'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_O'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Ne'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Mg'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Si'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_S'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Ca'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Fe'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
        except:
            ds.add_field(('star','metals'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
    else:
        ds.add_field(('star','metals'), function=_starmetals, sampling_type='particle',units="code_metallicity")

    if ('PartType0', 'Metallicity_00') in ds.derived_field_list:
        try:
            ds.add_field(('gas','metals'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_He'), function=_gasmetals_00, sampling_type='particle', units="code_metallicity")
            ds.add_field(('gas','metals_C'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_N'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_O'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Ne'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Mg'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Si'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_S'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Ca'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Fe'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
        except:
            ds.add_field(('gas','metals'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")

        ds.add_field(('metal','dens'), function=_metaldens_00, sampling_type='particle',units="g/cm**3")

        # we add this as part type 0 (a non-general name) as it gets
        # smoothed immediately and that's all we end up using downstream
        ds.add_field(('PartType0', 'metalmass'), function=_metalmass_00, sampling_type='particle',units="g")

    else:
        ds.add_field(('gas','metals'), function=_gasmetals, sampling_type='particle',units="code_metallicity")
        ds.add_field(('metal','dens'), function=_metaldens, sampling_type='particle',units="g/cm**3")
        ds.add_field(('PartType0', 'metalmass'), function=_metalmass, sampling_type='particle',units="g")

    if add_smoothed_quantities == True: ds.add_field(('metal','smoothedmasses'), function=_metalsmoothedmasses, sampling_type='particle',units='code_metallicity')

    ds.add_field(('gas','masses'), function=_gasmasses, sampling_type='particle',units='g')
    ds.add_field(('gas','fh2'), function=_gasfh2, sampling_type='particle',units='dimensionless')
    ds.add_field(('gas','sfr'), function=_gassfr, sampling_type='particle',units='g/s')
    ds.add_field(('gas','smoothinglength'),function=_gassmoothinglength,sampling_type='particle',units='pc')

    # get the dust mass
    ds.add_field(('dust','mass'), function=_dustmass_manual, sampling_type='particle',units='code_mass')
    ds.add_deposited_particle_field(("PartType0", "Dust_Masses"), "sum")
    #this just saves (redundantly) for passive dust 'manual' models the dust mass in 'particle_dust','mass' tuple.
    ds.add_field(('particle_dust','mass'),function=_return_dust_mass,units='code_mass',sampling_type='particle')
    
    if add_smoothed_quantities == True: ds.add_field(('dust','smoothedmasses'), function=_dustsmoothedmasses, sampling_type='particle',units='code_mass')
    
    ds.add_field(('star','masses'), function=_starmasses, sampling_type='particle',units='g')
    ds.add_field(('star','coordinates'), function=_starcoordinates, sampling_type='particle',units='cm')
    ds.add_field(('star','formationtime'), function=_starformationtime, sampling_type='particle',units='dimensionless')
    
    ds.add_field(('stellar','ages'),function=_stellarages,sampling_type='particle',units='Gyr')
    
    
    if add_smoothed_quantities == True: ds.add_field(('star','smoothedmasses'), function=_starsmoothedmasses, sampling_type='particle',units='g')

    ds.add_field(('gas','density'), function=_gasdensity, sampling_type='particle',units='g/cm**3')
    ds.add_field(('gas','coordinates'), function=_gascoordinates, sampling_type='particle',units='cm')
    if add_smoothed_quantities == True:
        ds.add_field(('gas','smootheddensity'), function=_gassmootheddensity, sampling_type='particle',units='g/cm**3')
        ds.add_field(('gas','smoothedmasses'), function=_gassmoothedmasses, sampling_type='particle',units='g')


        #try:
        ds.add_field(('gas','smoothedmetals'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_He'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_C'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_N'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_O'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Ne'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Mg'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Si'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_S'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Ca'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Fe'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        #except:
        #    ds.add_field(('gassmoothedmetals'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')


    return ds



def find_max_level(refined):
    #loop through refined and figure out how many trues we can have in a row
    master_max_level = 0
    max_level = 0.
    for i in range(len(refined)-1):
        if refined[i+1] == True:
            max_level+=1
        else:
            max_level = 0
        if max_level > master_max_level: master_max_level = max_level
        if max_level > 20: pdb.set_trace()
    return master_max_level



def octree_zoom_bbox_filter(fname,ds,bbox0):

    ds.index
    ad = ds.all_data()

    print ('\n\n')
    print ('----------------------------')
    print ("[octree zoom_bbox_filter:] Calculating Center of Mass")
    gas_com_x = np.sum(ad["gas","density"] * ad["gas","coordinates"][:,0])/np.sum(ad["gas","density"])
    gas_com_y = np.sum(ad["gas","density"] * ad["gas","coordinates"][:,1])/np.sum(ad["gas","density"])
    gas_com_z = np.sum(ad["gas","density"] * ad["gas","coordinates"][:,2])/np.sum(ad["gas","density"])
    com = [gas_com_x,gas_com_y,gas_com_z]

    print ("[octree zoom_bbox_filter:] Center of Mass is at coordinates (kpc): ",com)
    center = [x_cent,y_cent,z_cent]
    print ('[octree zoom_bbox_filter:] using center: ',center)

    box_len = 50
    box_len = ds.quan(box_len,'kpc')
    box_len = float(box_len.to('code_length').value)
    bbox_lim = box_len

    bbox1 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
    print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox1)

    ds = field_add(fname,bounding_box = bbox1,ds=ds,add_smoothed_quantities=True)
    reg = ds.region(center=center,left_edge = np.asarray(center)-bbox_lim,right_edge = np.asarray(center)+bbox_lim)
    
    left = np.array([pos[0] for pos in bbox1])
    right = np.array([pos[1] for pos in bbox1])
    octree = ds.octree(left, right, n_ref=nref)#, force_build=True)
    
    reg.parameters={}
    reg.parameters['octree'] = octree
    return reg


def gridstats(fc1,fw1):

    volume = (fw1**3)
    xmin = fc1[:,0]-fw1[:,0]/2.
    xmax = fc1[:,0]+fw1[:,0]/2.
    ymin = fc1[:,1]-fw1[:,1]/2.
    ymax = fc1[:,1]+fw1[:,1]/2.
    zmin = fc1[:,2]-fw1[:,2]/2.
    zmax = fc1[:,2]+fw1[:,2]/2.


    print ('----------------------------')
    print ('Grid Statistics')
    print ('----------------------------')

    print ('Smallest Cell Edge: ',np.min(fw1.in_units('kpc')))
    print ('Smallest Cell Volume: ',np.min(volume.in_units('kpc**3')))

    print ('Biggest Cell Edge: ',np.max(fw1.in_units('kpc')))
    print ('Biggest Cell Volume: ',np.max(volume.in_units('kpc**3')))

    print ('Left Edge: ',np.min(xmin.in_units('kpc')))
    print ('Right Edge: ',np.max(xmax.in_units('kpc')))
    print ('Bottom Edge: ',np.min(ymin.in_units('kpc')))
    print ('Top Edge: ',np.max(ymax.in_units('kpc')))
    print ('Nearest Edge: ',np.min(zmin.in_units('kpc')))
    print ('Farthest Edge: ',np.max(ymax.in_units('kpc')))


def yt_octree_generate(fname):

    print('[grid_construction]: bbox_lim = ', bbox_lim)

    bbox = [[-2.*bbox_lim, 2.*bbox_lim],
            [-2.*bbox_lim, 2.*bbox_lim],
            [-2.*bbox_lim, 2.*bbox_lim]]

    # load the DS and add pd fields.  this is the first field addition
    # of the simulation, so we don't yet need to add the smoothed
    # quantities (which can take some time in yt4.x).
    ds = field_add(fname, bounding_box=bbox,add_smoothed_quantities=False)

    #now zoom in.
    reg = octree_zoom_bbox_filter(fname, ds, bbox)


    refined = reg.parameters['octree'][('index','refined')].astype('bool')
     
    xpos = (reg.parameters['octree'][('index','x')])[~refined]
    ypos = (reg.parameters['octree'][('index','y')])[~refined]
    zpos = (reg.parameters['octree'][('index','z')])[~refined]
    
    #comebine these into the fc1 array with dimensions (nparticles,3)
    fc1 = np.array([xpos,ypos,zpos]).T
    
    dx = reg.parameters['octree'][('index','dx')][~refined]
    dy = reg.parameters['octree'][('index','dy')][~refined]
    dz = reg.parameters['octree'][('index','dz')][~refined]
    fw1 = np.array([dx,dy,dz]).T 
    
    n_ref = reg.parameters['octree'].n_ref
    #max_level = find_max_level(refined)#'max_level not yet implemented in yt4.x'
    #note, we could figure this out from the max number of trues in a row
    nocts = len(refined)-np.sum(refined)
    # convert fc1 and fw1 to YTArrays
    fc1 = ds.arr(fc1, 'code_length')
    fw1 = ds.arr(fw1, 'code_length')


    print('----------------------------')
    print('yt Octree Construction Stats')
    print('----------------------------')
    print(' n_ref = ', n_ref)
    #print(' max_level = ', max_level)
    print(' nocts = ', nocts)
    print('----------------------------')

    gridstats(fc1, fw1)

    refinements = 2**(3*0)
    refined2 = []
    for r in refined:
        if r == 1:
            refined2.append(True)
        if r == 0:
            refined2.append(np.zeros(refinements).astype('bool'))
    refined = np.hstack(refined2)


    dust_smoothed_manual = manual_oct(reg, refined)
    dust_smoothed = dust_smoothed_manual


    return refined, dust_smoothed, fc1, fw1, reg, ds



def hyperion_octree_stats(refined):

    if not (len(refined) - 1) % 8 == 0:
        raise ValueError("refined should have shape 8 * n + 1")

    refined = np.array(refined, dtype=bool, copy=False)

    print("Number of items in refined            : {0}".format(len(refined)))
    print("Number of True values in refined      : {0}".format(np.sum(refined)))
    print("Number of False values in refined     : {0}".format(np.sum(~refined)))

    def check_recursive(refined, current_i=0, max_level=0):

        if refined[current_i]:
            current_i += 1
            max_levels = []
            for i in range(8):
                current_i, max_level_indiv = check_recursive(refined, current_i, max_level+1)
                max_levels.append(max_level_indiv)
            max_level = max(max_levels)
        else:
            current_i += 1
        return current_i, max_level

    try:
        final_i, max_level = check_recursive(refined)
    except IndexError as e:
        max_level = 'unknown'
        consistent = False
    else:
        consistent = True

    print("Array is self-consistent for Hyperion : {0}".format("yes" if consistent else "no"))
    print("Maximum number of levels              : {0}".format(max_level))

    return max_level

def sanity_check(refined_string,max_level):
    print ('Entering the Octree Sanity Check')
    content = refined_string
    prev = content

    print ('inside sanity_check: max_level = ',max_level)

    for iter in range(max_level):

        content = content.replace('TFFFFFFFF', 'F')

        print ("Length: {0}".format(len(content)))
        if content == prev:
            break

            prev = content


    print ('Converged: {0} (should be F)'.format(content))


def test_octree(refined,max_level):
    #convert refined to a string to make life easy for octree_sanity_check
    rs = str(refined)

    rs = rs.replace("[","").replace(',','').replace("]","")
    rs = rs.replace("True","T").replace("False","F")
    rs = rs.replace(" ","")



    sanity_check(rs,max_level)


def dump_cell_info(refined,fc1,fw1,xmin,xmax,ymin,ymax,zmin,zmax):
    outfile = f"/orange/narayanan/d.zimmerman/camels_scripts/cell_info.33_{galaxy_num}.npz"
    np.savez(outfile,refined=refined,fc1=fc1,fw1=fw1,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax)

def find_order(refined):
    """
    Find the index array to use to sort the ``refined`` and ``density`` arrays.
    """

    order = np.zeros(refined.shape)

    if not refined[0]:
        return [0]


    def find_nested(i):
        cells = [i]
        for cell in range(8):
            i += 1

            if i >= len(refined):
                return i,np.hstack(cells)

            if refined[i]:
                parent = i
                i, sub_cells = find_nested(i)
                cells.append(sub_cells)
            else:
                cells.append(i)
        cells = [cells[j] for j in [0,1,5,3,7,2,6,4,8]]
        return i, np.hstack(cells)

    return find_nested(0)[1]

def sph_m_gen(fname):

    refined,dustdens,fc1,fw1,reg,ds = yt_octree_generate(fname)

    xmin = (fc1[:,0]-fw1[:,0]/2.).to('cm') #in proper cm
    xmax = (fc1[:,0]+fw1[:,0]/2.).to('cm')
    ymin = (fc1[:,1]-fw1[:,1]/2.).to('cm')
    ymax = (fc1[:,1]+fw1[:,1]/2.).to('cm')
    zmin = (fc1[:,2]-fw1[:,2]/2.).to('cm')
    zmax = (fc1[:,2]+fw1[:,2]/2.).to('cm')

    #dx,dy,dz are the edges of the parent grid
    dx = (np.max(xmax)-np.min(xmin)).value
    dy = (np.max(ymax)-np.min(ymin)).value
    dz = (np.max(zmax)-np.min(zmin)).value


    xcent = float(ds.quan(x_cent,"code_length").to('cm').value)
    ycent = float(ds.quan(y_cent,"code_length").to('cm').value)
    zcent = float(ds.quan(z_cent,"code_length").to('cm').value)

    boost = np.array([xcent,ycent,zcent])
    print ('[sph_tributary] boost = ',boost)
    print ('[sph_tributary] xmin (pc)= ',np.min(xmin.to('pc')))
    print ('[sph_tributary] xmax (pc)= ',np.max(xmax.to('pc')))
    print ('[sph_tributary] ymin (pc)= ',np.min(ymin.to('pc')))
    print ('[sph_tributary] ymax (pc)= ',np.max(ymax.to('pc')))
    print ('[sph_tributary] zmin (pc)= ',np.min(zmin.to('pc')))
    print ('[sph_tributary] zmax (pc)= ',np.max(zmax.to('pc')))
    #Tom Robitaille's conversion from z-first ordering (yt's default) to
    #x-first ordering (the script should work both ways)

    refined_array = np.array(refined)
    refined_array = np.squeeze(refined_array)

    order = find_order(refined_array)
    refined_reordered = []
    dustdens_reordered = np.zeros(len(order))

    for i in range(len(order)):
        refined_reordered.append(refined[order[i]])
        dustdens_reordered[i] = dustdens[order[i]]

    refined = refined_reordered
    dustdens=dustdens_reordered

    #hyperion octree stats
    max_level = hyperion_octree_stats(refined)


    test_octree(refined,max_level)


    #save some information that can be used in the PAH model compute
    #an effective 'size' of a cell by density = mass/volume and assume
    #spherical geometry.  similarly, saving the particle location information
    dump_cell_info(refined,fc1.to('cm'),fw1.to('cm'),xmin,xmax,ymin,ymax,zmin,zmax)
    reg.parameters['cell_size']=fw1.convert_to_units('cm') #so that we can have a uniform naming scheme for different front ends for saving in analytics/dump_data(
    reg.parameters['cell_position'] = fc1

    np.save('refined.npy',refined)
    np.save('density.npy',dustdens)


    #save some information that can be used in the PAH model
    reg.parameters['fw1'] = fw1.convert_to_units('cm')


    #========================================================================
    #Initialize Hyperion Model
    #========================================================================

    m = Model()

    #save in the m__dict__ that we're in an oct geometry
    #m.__dict__['grid_type']='oct'

    print ('Setting Octree Grid with Parameters: ')



    #m.set_octree_grid(xcent,ycent,zcent,
    #                  dx,dy,dz,refined)
    #m.set_octree_grid(0,0,0,dx/2,dy/2,dz/2,refined)


    #energy_density_absorbed=energy_density_absorbed_by_CMB()
    #specific_energy = np.repeat(energy_density_absorbed.value,dustdens.shape)


    #d = SphericalDust(dustdir+dustfile)
    #m.add_density_grid(dustdens,d,specific_energy=specific_energy)
    #m.set_specific_energy_type('additional')

    return m,xcent,ycent,zcent,dx,dy,dz,reg,ds,boost

def dump_data(reg,model):

    particle_fh2 = reg["gas","fh2"]
    particle_fh1 = np.ones(len(particle_fh2))-particle_fh2
    particle_gas_mass = reg["gas","masses"]
    particle_star_mass = reg["star","masses"]
    particle_star_metallicity = reg["star","metals"]
    #particle_stellar_formation_time = reg["starformationtime"]
    particle_stellar_formation_time = reg["stellar","ages"]
    particle_sfr = reg['gas','sfr'].in_units('Msun/yr')


    try: cell_size = reg.parameters['cell_size']
    except: cell_size = -1

    #save dustmasses.  for particle type codes where the particles are
    #projected onto an octree (gizmo/gadget) this will be the
    #smoothedmasses field; else, it will just be the ('dust','mass')
    #tuple which is the mesh values
    try: grid_dustmass = reg['dust','smoothedmasses'].in_units('Msun')
    except: grid_dustmass = reg["dust","mass"].in_units('Msun')  #for arepo/AMR codes,
                                                 #this is the grid
                                                 #values (and also
                                                 #parttype0 for
                                                 #arepo).  for
                                                 #gizmo/gadget, this
                                                 #is particle
                                                 #information which is
                                                 #why it has to be in
                                                 #the except
                                                 #statement, and not
                                                 #the try statemetn in logical flow.

    #if we have separate particle dust masses, save them.
    try: particle_dustmass = reg['particle_dust','mass'].in_units('Msun')
    except: particle_dustmass = -1


    #these are in try/excepts in case we're not dealing with gadget and yt 3.x

    try:
        grid_gas_mass = reg["gas","smoothedmasses"]
    except: grid_gas_mass = -1
    #for sph we have a different nomenclature thanks to the smoothing onto an octree
    try:
        grid_gas_metallicity = []
        grid_gas_metallicity.append(reg["gas","smoothedmetals"].value)
        abund_el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
        for i in abund_el:
            grid_gas_metallicity.append(reg["gas","smoothedmetals_"+str(i)].value)

    except: grid_gas_metallicity = reg["gas","smoothedmetals"].value

    try: grid_star_mass = reg["star","smoothedmasses"]
    except: grid_star_mass = -1

    #PAH information
    try: grid_PAH_luminosity = reg.parameters['grid_PAH_luminosity']
    except: grid_PAH_luminosity = -1
    try: grid_neutral_PAH_luminosity = reg.parameters['grid_neutral_PAH_luminosity']
    except: grid_neutral_PAH_luminosity = -1
    try: grid_ion_PAH_luminosity = reg.parameters['grid_ion_PAH_luminosity']
    except: grid_ion_PAH_luminosity = -1


    try: PAH_lam = reg.parameters['PAH_lam']
    except: PAH_lam = -1
    try: total_PAH_luminosity = reg.parameters['total_PAH_luminosity']
    except: total_PAH_luminosity = -1
    try: total_neutral_PAH_luminosity = reg.parameters['total_neutral_PAH_luminosity']
    except: total_neutral_PAH_luminosity = -1
    try: total_ion_PAH_luminosity = reg.parameters['total_ion_PAH_luminosity']
    except: total_ion_PAH_luminosity = -1


    try: integrated_grid_PAH_luminosity = reg.parameters['integrated_grid_PAH_luminosity']
    except: integrated_grid_PAH_luminosity = -1
    try: integrated_grid_neutral_PAH_luminosity = reg.parameters['integrated_grid_neutral_PAH_luminosity']
    except: integrated_grid_neutral_PAH_luminosity = -1
    try: integrated_grid_ion_PAH_luminosity = reg.parameters['integrated_grid_ion_PAH_luminosity']
    except: integrated_grid_ion_PAH_luminosity = -1



    try: q_pah = reg.parameters['q_pah']
    except: q_pah = -1
    try: particle_mass_weighted_gsd = reg.parameters['particle_mass_weighted_gsd']
    except: particle_mass_weighted_gsd = -1
    try: grid_mass_weighted_gsd = reg.parameters['grid_mass_weighted_gsd']
    except: grid_mass_weighted_gsd = -1
    try: simulation_sizes = reg.parameters['simulation_sizes']
    except: simulation_sizes = -1
    try: beta_nnls = reg.parameters['beta_nnls']
    except: beta_nnls = -1
    #get tdust
    #m = ModelOutput(model.outputfile+'.sed')
    #oct = m.get_quantities()
    #tdust_ds = oct.to_yt()
    #tdust_ad = tdust_ds.all_data()
    #tdust = tdust_ad[ ('gas', 'temperature')]


    #try:
    outfile = f"/orange/narayanan/d.zimmerman/camels_scripts/grid_physical_properties.33_galaxy{galaxy_num}.npz"
    #except:
     #   outfile = cfg.model.PD_output_dir+"/grid_physical_properties."+cfg.model.snapnum_str+".npz"

    np.savez(outfile,particle_fh2=particle_fh2,particle_fh1 = particle_fh1,particle_gas_mass = particle_gas_mass,particle_star_mass = particle_star_mass,particle_star_metallicity = particle_star_metallicity,particle_stellar_formation_time = particle_stellar_formation_time,grid_gas_metallicity = grid_gas_metallicity,grid_gas_mass = grid_gas_mass,grid_star_mass = grid_star_mass,particle_sfr = particle_sfr,particle_dustmass = particle_dustmass,grid_dustmass=grid_dustmass,grid_PAH_luminosity = grid_PAH_luminosity,grid_neutral_PAH_luminosity=grid_neutral_PAH_luminosity,grid_ion_PAH_luminosity=grid_ion_PAH_luminosity,PAH_lam=PAH_lam,total_PAH_luminosity = total_PAH_luminosity,total_neutral_PAH_luminosity=total_neutral_PAH_luminosity,total_ion_PAH_luminosity=total_ion_PAH_luminosity,integrated_grid_PAH_luminosity = integrated_grid_PAH_luminosity,integrated_grid_neutral_PAH_luminosity=integrated_grid_neutral_PAH_luminosity,integrated_grid_ion_PAH_luminosity=integrated_grid_ion_PAH_luminosity,q_pah=q_pah,particle_mass_weighted_gsd = particle_mass_weighted_gsd,grid_mass_weighted_gsd = grid_mass_weighted_gsd,simulation_sizes=simulation_sizes,cell_size=cell_size,beta_nnls=beta_nnls)#,tdust = tdust)




def manual_oct(reg,refined):
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]


    density_smoothed = reg["gas","smootheddensity"]
    metallicity_smoothed = reg["gas","smoothedmetals"]
    masses_smoothed = reg["gas","smoothedmasses"]


    try:
        smoothed_dust_masses = reg[('dust','smoothedmasses')]
    except:
        raise KeyError('Dust mass information not present in this snapshot. Please set another dust grid type in the parameters.')
    dust_to_gas_ratio = smoothed_dust_masses.in_units('g')/masses_smoothed
    #masses_smoothed can be 0 at some places; this will make dtg nan
    #out even though it would eventually get multiplied to 0 when we
    #multiply by density smoothed.  So, to be stable we nan_to_num it.
    dust_to_gas_ratio = np.nan_to_num(dust_to_gas_ratio)
    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed
    return dust_smoothed








bbox = [[-2.*bbox_lim,2.*bbox_lim],
         [-2.*bbox_lim,2.*bbox_lim],
         [-2.*bbox_lim,2.*bbox_lim]]
ds = yt.load(fname,bounding_box = bbox)
#ds.index
#print ('[front_end_controller:] bounding_box being used')


#ds = field_add(fname,bounding_box=bbox)
#m,xc,yc,zc,dx,dy,sz,reg,ds2,boost = sph_m_gen(fname)
#dump_data(reg, model)



center = [x_cent,y_cent,z_cent]
print ('[octree zoom_bbox_filter:] using center: ',center)
box_len = 50
box_len = ds.quan(box_len,'kpc')
box_len = float(box_len.to('code_length').value)
bbox_lim = box_len
bbox2 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
print()
print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox2)

ds2 = yt.load(fname,bounding_box = bbox2)
ds2.index

left = np.array([pos[0] for pos in bbox2])
right = np.array([pos[1] for pos in bbox2])
#octree = ds.octree(left, right, over_refine_factor=cfg.par.oref, n_ref=cfg.par.n_ref, force_build=True)
octree = ds2.octree(left,right,n_ref=nref)
ds2.parameters['octree'] = octree
reg = ds2.all_data()

print(octree)
#octree.print_all_nodes()
max_rho_xyz = reg.argmax(("gas", "density"))
print(max_rho_xyz)
print(octree["PartType0","Masses"].to('Msun'))
print(np.sum(octree["PartType0","Masses"].to('Msun')))
val0 = octree["PartType0","Density"].to('Msun/kpc**3')
refined = octree[('index','refined')].astype('bool')
volume = octree['index','dx'][~refined]*octree['index','dy'][~refined]*octree['index','dz'][~refined]
print(np.sum(volume.to('kpc**3')*val0))
