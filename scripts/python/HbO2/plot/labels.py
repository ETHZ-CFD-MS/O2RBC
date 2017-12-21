"""Function to set axis labels for Hb-O2 simulation."""

import warnings

import matplotlib.pyplot as plt

mmHg = r'\mathrm{mm\,Hg}'
mlO2 = r'\mathrm{ml\,O_2}'
umO2 = r'\mathrm{\muup m^3\,O_2}'
IVR = mmHg + r'\mathrm{\:\muup m\,s/}' + umO2
m0_unit = r'\mathrm{10^{-3}\muup m^3\,O_2\,\muup m^{-3}\,s^{-1}}'
jt_unit = r'\mathrm{\muup m^3\,O_2\,\muup m^{-1}\,s^{-1}}'

OFVarToLatex = {'Hb_mean':           r'S',
                'deltaLD':           r'\mathrm{LD}',
                'deltaS':            r'\Delta S_a',
                'deltaV':            r'\Delta v_{\mathrm{rbc}}',
                'H':                 r'\mathrm{capillary\;spacing}',
                'LD':                r'\mathrm{LD}',
                'LDMean':            r'\mathrm{LD}',
                'O2ConsumptionRate': r'M_0',
                'radiusRBC':         r'r_{\mathrm{rbc}}',
                'radiusPlasma':      r'r_p',
                'sigmaS':            r'\sigma_{S,a}',
                'RBCVelocity':       r'v_{\mathrm{rbc}}'}

OFVarToUnit = {'Hb_mean':           '',
               'deltaS':            '',
               'deltaV':            'mm/s',
               'H':                 'um',
               'LD':                '',
               'LDMean':            '',
               'O2ConsumptionRate': m0_unit,
               'radiusRBC':         'um',
               'radiusPlasma':      'um',
               'sigmaS':            '',
               'RBCVelocity':       'mm/s',
               }

OFVarScaling = {'deltaV':            1e3,
                'H':                 1e6,
                'O2ConsumptionRate': 1e-3,
                'radiusRBC':         1e6,
                'radiusPlasma':      1e6,
                'RBCVelocity':       1e3}

varToLatex = {'LD':   r'\mathrm{LD}',
              'PO2':  r'\mathrm{PO}_2',
              'vrbc': r'v_{\mathrm{rbc}}',
              'MTC':  r'\mathrm{MTC}'}

unitToLatex = {'mmHg':   mmHg,
               'um':     r'\mathrm{\muup m}',
               's':      r'\mathrm{s}',
               'ms':     r'\mathrm{ms}',
               'us':     r'\mathrm{\muup s}',
               'mm/s':   r'\mathrm{mm\,s^{-1}}',
               'M_0':    m0_unit,
               'IVR':    IVR,  # unit of intravascular resistance
               '1e6IVR': r'10^6\,' + IVR,
               'MTC':    r'10^{-6}\mathrm{ml\,O2\,s^{-1}\,mm\,Hg^{-1}\,cm^{-2}}',
               'j_t':    jt_unit}


def varUnitString(var, unit):
    """Construct a formatted string with the given variable and unit.

    Both arguments are looked up in the dictionaries 'varToLatex' and 
    'unitToLatex' for argument conversion. If they are not found, these
    arguments will be used as given. If unit is None, no unit will be used
    and the square brackets will be dropped.

    Return: 
        a string of the form '<var> [<unit>]' or '<var>'
    """
    try:
        latexVar  = varToLatex[var] 
    except KeyError:
        latexVar = var
    try:
        latexUnit = unitToLatex[unit]
    except KeyError:
        latexUnit = unit
    s = r'%s' % latexVar
    if unit:
        s += r' \; [%s]' % latexUnit
    return r'$%s$' % s


def setTitle(var, unit):
    plt.title(varUnitString(var, unit))


def setLabel(axis, var, unit, **kwargs):
    if axis not in ['x', 'y']:
        raise ValueError('Invalid axis %s' % str(axis))
    if 'ax' in kwargs:
        ax = kwargs.get('ax')
        labelFunction = getattr(ax, 'set_%slabel' % axis)
    else:
        labelFunction = getattr(plt, '%slabel' % axis)
    labelFunction(varUnitString(var, unit))


def setXLabel(var, unit=None, **kwargs):
    setLabel('x', var, unit, **kwargs)


def setYLabel(var, unit=None, **kwargs):
    setLabel('y', var, unit, **kwargs)


def openFOAMVarNameToLatex(varName):
    """Return the LaTeX version of a variable name used in OpenFOAM."""
    try:
        return OFVarToLatex[varName]
    except KeyError:
        warnings.warn("""Variable name '{}' not recognized, 
                         returning empty string""".format(varName))
        return ""


def openFOAMVarNameToUnit(varName):
    """Return the unit of a variable name used in OpenFOAM."""
    try:
        return OFVarToUnit[varName]
    except KeyError:
        warnings.warn("""Variable name '{}' not recognized, 
                         returning empty string""".format(varName))
        return ""


def openFOAMVarScaling(varName):
    """Return the plotting scaling factor for a variable name used in OpenFOAM."""
    try:
        return OFVarScaling[varName]
    except KeyError:
        return 1
