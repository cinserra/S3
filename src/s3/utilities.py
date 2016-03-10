####################################
def cardelli(lam, Av, Rv=3.1, Alambda = True, debug=False):
    import numpy as np
    """ 
    keywords:
        Alambda        <bool>  returns +2.5*1./log(10.)*tau
        Av        <float>    extinction value
        Rv        <float> extinction param. (def: 3.1)
    """
    if type(lam) == float:
        _lam = np.asarray([lam])
    else:
        _lam = lam[:]

    #initialisation of variables
    x = 10000./(_lam) #to convert from micron to Angstrom
    a = np.zeros(np.size(x))
    b = np.zeros(np.size(x))
    #Infrared (Eq 2a,2b)
    inst = np.where ((x >= 0.3) & (x < 1.1))
    a[inst] =  0.574*x[inst]**1.61
    b[inst] = -0.527*x[inst]**1.61
    ##############  Optical & Near IR ###############
    #Eq 3a, 3b
    inst = np.where ((x >= 1.1) & (x <= 3.3))
    y = x[inst]-1.82
    a[inst] = 1. + 0.17699*y   - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7
    b[inst] =      1.41338*y   + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7
    ############# UV  ############################### 
    #Eq 4a, 4b
    inst = np.where ((x >= 3.3) & (x <= 8.0))
    a[inst] =  1.752 - 0.316*x[inst] - 0.104/((x[inst]-4.67)**2+0.341)
    b[inst] = -3.090 + 1.825*x[inst] + 1.206/((x[inst]-4.62)**2+0.263)

    inst = np.where ((x >= 5.9) & (x <= 8.0))
    Fa     = -0.04473*(x[inst]-5.9)**2 - 0.009779*(x[inst]-5.9)**3
    Fb     =  0.21300*(x[inst]-5.9)**2 + 0.120700*(x[inst]-5.9)**3
    a[inst] = a[inst] + Fa
    b[inst] = b[inst] + Fb
    ########## Far UV ##############################
    #Eq 5a, 5b
    inst = np.where ((x >= 8.0) & (x <= 10.0))
    #Fa = Fb = 0
    a[inst] = -1.073 - 0.628*(x[inst]-8.) + 0.137*((x[inst]-8.)**2) - 0.070*(x[inst]-8.)**3 
    b[inst] = 13.670 + 4.257*(x[inst]-8.) + 0.420*((x[inst]-8.)**2) + 0.374*(x[inst]-8.)**3

    # Case of -values x out of range [0.3,10.0]
    inst = np.where ((x > 10.0) | (x < 0.3))
    a[inst] = 0.0
    b[inst] = 0.0

    #Return Extinction vector 
    #Eq 1
    if (Alambda == True):
        return 10**(0.4*( 2.5*1./np.log(10.)*( a + b/Rv ) * Av)) #added the factor 10^(0.4*cardelli result) to work with spectra
    else:
        return 10**(0.4*( ( a + b/Rv ) * Av)) #added the factor 10^(0.4*cardelli result) to work with spectra

##############################################
def bbody(lam,T,A):
    import numpy as np
    """ 
    keywords:
        no keywords
    """
    Blam = A*(2*6.6261e-27*(2.9979e10)**2/(lam*1e-8)**5)/(np.exp(6.6261e-27*2.9979e10/(lam*1e-8*1.3806e-16*T))-1)
    return Blam

# units:
# Blam in erg/(s cm^3)
# speed of light c = 2.9979e10 cm/s
# Planck constant h = 6.6261e-27 erg*s
# wavelengths lam input in Angstrom then converted automatically to cm
# Boltzmann constant kB = 1.3806e-16 erg/K
# T in K

##############################################
def conv(wx,fy,wxf,tran,wspm,wspmx,wfm,wfmx):
    from scipy import *
    import numpy as np
    from scipy.interpolate import interp1d
    if len(wx) != len(wxf):
        if len(wx) > len(wxf):
            f = interp1d(wxf,tran,bounds_error=False,fill_value=0.)
            wxf = np.linspace(wspm,wspmx,len(wx))
            tran = f(wxf)
        if len(wxf) > len(wx):
            f = interp1d(wx,fy,bounds_error=False,fill_value=0.)
            wx = np.linspace(wfm,wfmx,len(wxf))
            fy = f(wx)
    cf = fy * tran
    return cf, wx #return flux and wavelength for the integration after the convolution
# interpolating the two responses to match the length and sampling coverage
##############################################