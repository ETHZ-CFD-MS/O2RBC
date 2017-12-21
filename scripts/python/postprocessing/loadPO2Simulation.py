import ConfigParser
import os
import STLUtils


# Read the parameter file multi_RBC_params and return data in a
# dictionary.
def loadMultiRBCParams(caseName):
    param_file = "multi_RBC_params.txt"
    config = ConfigParser.RawConfigParser()
    config.read(param_file)

    # extract data from config file
    nRBC             = config.getint('Main', 'nRBC')
    line_density_str = config.get('Main', 'line_density')
    Rt               = config.getfloat('Main', 'Rt')
    STLFile          = config.get('Main', 'STLFile')
    PO2_RBC_str      = config.get('Main', 'PO2_RBC')
    PO2_tissue       = config.getfloat('Main', 'PO2_tissue')
    monitored_RBC    = config.getint('Main', 'monitored_RBC')
    dx               = config.getfloat('Main', 'dx')
    dy               = config.getfloat('Main', 'dy')

    # split the string with PO2_RBC
    PO2_RBC = [float(x) for x in PO2_RBC_str.split(',')]
    if len(PO2_RBC) != nRBC:
        print "Error, the number of PO2_RBC values does not match " \
            + "nRBC."
        sys.exit()

    # split the string with the line density/ies
    line_density = [float(x) for x in line_density_str.split(',')]
    if len(line_density) != nRBC and len(line_density) != 1:
        print "Error, the number of values of line density is not equal to " \
            + "1 or nRBC."
        sys.exit()

    # if a constant line density is given, fill a list of size nRBC with
    # this value.
    if (len(line_density) == 1):
        line_density = [line_density[0] for i in xrange(nRBC)]

    # extract RBC length
    if os.path.isfile(STLFile):
        L_RBC = STLUtils.getLength(STLFile, 0)
    else:
        L_RBC = 8.3467925710416212e-06 # length for RBC cylinder, V_RBC = 59e-18, r_RBC = 1.5e-6

    # return all data stored in a dictionary
    return dict([('nRBC', nRBC), \
                 ('line_density', line_density), \
                 ('Rt', Rt), \
                 ('STLFile', STLFile), \
                 ('PO2_RBC', PO2_RBC), \
                 ('PO2_tissue', PO2_tissue), \
                 ('monitored_RBC', monitored_RBC), \
                 ('dx', dx), \
                 ('dy', dy), \
                 ('L_RBC', L_RBC)])

# Read the parameter file sim_params and return data in a
# dictionary.
def loadSimParams(caseName):
    param_file = "sim_params.txt"
    config = ConfigParser.RawConfigParser()
    config.read('%s/%s' % (caseName, param_file))

    # extract data from config file
    nRBC            = config.getint('Main', 'nRBC')
    LDomain         = config.getfloat('Main', 'LDomain')
    lineDensity     = config.getfloat('Main', 'lineDensity')
    CFL             = config.getfloat('Main', 'CFL')
    RtLeft          = config.getfloat('Main', 'RtLeft')
    RtRight         = config.getfloat('Main', 'RtRight')
    Rp              = config.getfloat('Main', 'Rp')
    Rw              = config.getfloat('Main', 'Rw')
    PO2RBCInlet     = config.getfloat('Main', 'PO2RBCInlet')
    dx              = config.getfloat('Main', 'dx')
    dy              = config.getfloat('Main', 'dy')
    dxLag           = config.getfloat('Main', 'dxLag')
    dyLag           = config.getfloat('Main', 'dyLag')
    nyTissue        = config.getint('Main', 'nyTissue')

    # extract RBC length
    L_RBC = 8.3467925710416212e-06 # length for RBC cylinder, V_RBC = 59e-18, r_RBC = 1.5e-6
    L_RBC_centerline = L_RBC
    STLFile = ""
    try:
        STLFile     = config.get('Main', 'STLFile')
        if os.path.isfile(STLFile):
            STL_path = "%s/%s" % (caseName, STLFile)
            L_RBC            = STLUtils.getLength(STL_path, 0)
            L_RBC_centerline = STLUtils.getCenterlineLength(STL_path)
    except ConfigParser.NoOptionError:
        pass

    # return all data stored in a dictionary
    return dict([('nRBC', nRBC), \
                 ('LDomain', LDomain), \
                 ('lineDensity', lineDensity), \
                 ('CFL', CFL), \
                 ('RtLeft', RtLeft), \
                 ('RtRight', RtRight), \
                 ('Rp', Rp), \
                 ('Rw', Rw), \
                 ('STLFile', STLFile), \
                 ('PO2RBCInlet', PO2RBCInlet), \
                 ('dx', dx), \
                 ('dy', dy), \
                 ('dxLag', dxLag), \
                 ('dyLag', dyLag), \
                 ('nyTissue', nyTissue), \
                 ('L_RBC', L_RBC), \
                 ('L_RBC_centerline', L_RBC_centerline)])

# Return a list which contains a list of times with OpenFOAM output.
def returnOutputTimes(caseName):
    times = []
    # extract time list for one simulation
    for f in os.listdir(caseName):
        if f.startswith('0') and f != '0.org':
            times.append(float(f))

    times.sort()
    return times

