#!/usr/bin/env python

#everything is needed to perform the script and maybe something else
from numpy import *
from scipy import *
from scipy import integrate
from scipy.interpolate import interp1d
import pyfits
import os,sys,string,shutil,math
from pylab import *
from scipy.optimize import curve_fit
import s3 #import metadata 
from s3.utilities import *  #import definitions
# pre-set plot parameters, resolution untouched since it is not needed (default=80 dpi) 
from pylab import rcParams
rcParams['figure.figsize'] = 11, 8
rcParams['figure.subplot.top'] = 0.95
rcParams['figure.subplot.right'] = 0.95
rcParams['figure.subplot.left'] = 0.11
###########################################
pypath = os.path.expandvars('$HOME')           # it copies login.cl if it is not in the same directory
if not os.path.isfile('login.cl'):
    shutil.copyfile(pypath+'/iraf/login.cl','login.cl')
###########################################

################### for the help ##################
from optparse import OptionParser

description = " K-correction script for a single, flux calibrated spectrum "
usage = "%prog "
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog " + str(s3.__version__))
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    parser.add_option("-r", "--redshifterr",dest="redshifterr", action="store", type="float" ,default=None,
                  help='Change the default error on your redshift (+/- 0.005) to estimate the K-correction errors')
    option,args = parser.parse_args()

###### moved here because OptionParser --version conflicts with pyraf version########
#what we need from iraf
from pyraf import iraf

########### Options that can be changed thanks to option parser #########
if option.redshifterr == None:
    _redshifterr = 0.005
else:
    _redshifterr = option.redshifterr
################ internal description #############

h="######################################################################\n"+\
  "#########  SuperNova Algorithm for K-correction Evaluation  ##########\n"+\
  "##################           S.N.A.K.E.           ####################\n"+\
  "##########          C. Inserra  v1.1.0 29/10/2015          ###########\n"+\
  "######################################################################\n"+\
  " K-correction based on the formula m(x) = M(y) + DM + K(y,x)\n"+\
  " BE SURE that the spectrum is flux calibrated  \n"+ \
  " If you use this code and find it useful, please give a thought \n"+ \
  " to cite it. \n"+ \
  " The reference is Inserra et al. 2015, ApJ submitted \n"+\
  "######################################################################\n"

print h 

#the path where the metatabs data are
filterdir=s3.__path__[0]+'/metadata/' # To set the directory where are the synphot tabs created

# cleaning process
os.system('rm -rf sn.txt')
os.system('rm -rf sn.fits')
os.system('rm -rf sn_xbbody.txt')
os.system('rm -rf sn_dez_xbbody.txt')
os.system('rm -rf sn_dez_xbbody_err1.txt')
os.system('rm -rf sn_dez_xbbody_err2.txt')
os.system('rm -rf bbody_sn_dez_fit.dat')
os.system('rm -rf bbody_sn_dez_fit_err1.dat')
os.system('rm -rf bbody_sn_dez_fit_err2.dat')
os.system('rm -rf bbody_sn_dez_fit.fits')
os.system('rm -rf bbody_sn_dez_fit_err1.fits')
os.system('rm -rf bbody_sn_dez_fit_err2.fits')
os.system('rm -rf bbody_sn_fit.fits')
os.system('rm -rf bbody_sn_fit.dat')
os.system('rm -rf sn_dez.txt')
os.system('rm -rf sn_dez_err1.txt')
os.system('rm -rf sn_dez_err2.txt')
os.system('rm -rf sn_dez.fits')
os.system('rm -rf sn_dez_err1.fits')
os.system('rm -rf sn_dez_err2.fits')
os.system('rm -rf sn_dez_dered.fits')
os.system('rm -rf sn_dez_dered_err1.fits')
os.system('rm -rf sn_dez_dered_err2.fits')
os.system('rm -rf sn_dez_dered.txt')
os.system('rm -rf sn_dez_dered_err1.txt')
os.system('rm -rf sn_dez_dered_err2.txt')
os.system('rm -rf sn_galdered_dez.fits')
os.system('rm -rf sn_galdered_dez_err1.fits')
os.system('rm -rf sn_galdered_dez_err2.fits')
os.system('rm -rf sn_galdered.fits')
os.system('rm -rf bsn_combo_dez.fits')
os.system('rm -rf bsn_combo_dez_err1.fits')
os.system('rm -rf bsn_combo_dez_err2.fits')
os.system('rm -rf bsn_combo_dez.txt')
os.system('rm -rf bsn_combo_dez_err1.txt')
os.system('rm -rf bsn_combo_dez_err2.txt')
os.system('rm -rf bsn_combo.fits')
os.system('rm -rf bsn_combo.txt')


#######################################################
# Variable definitions
#######################################################

print ''
print '#################################'
print '#   Available filters and ID    #'
print '#-------------------------------#'
print '# BESSEL:  U  B  V  R  I        #'
print '# SLOAN:   us gs rs is zs       #'
print '# UV:  NUV  FUV  uw2  um2  uw1  #'
print '# NIR:     J  H  K              #' 
print '#################################'
print ''  


filter1 = raw_input('Which is the observed filter FUV,NUV,uvw2,uvm2,uvw1,U,B,V,R,I,us,gs,rs,is,zs,J,H,K [rs] ? ')
if not filter1:
    filter1 = 'rs'

print ''

if filter1 == 'K' or filter1 == 'uvw1' or filter1 == 'uvw2' or filter1 == 'uvm2' or filter1 == 'H' or filter1 =='J':
    questionsys1 = raw_input('Do you want to use the Vega or the AB system ([vega],ab) ? ')
    if not questionsys1:
        questionsys1 = 'vega'

    if questionsys1 == 'ab': #assign the new file with the ZP in AB magnitude
        if filter1 == 'K':
            filter1 = 'K_ab'
        elif filter1 == 'J':
            filter1 = 'J_ab'
        elif filter1 == 'H':
            filter1 = 'H_ab'
        elif filter1 == 'uvw1':
            filter1 = 'uvw1_ab'
        elif filter1 == 'uvw2':
            filter1 = 'uvw2_ab'
        elif filter1 == 'uvm2':
            filter1 ='uvm2_ab'

print ''
_redshift = raw_input('What is the redshift [0.1] ? ')
if not _redshift:
    redshift = 0.1
else: 
    redshift = float(_redshift)

print ''

#######################################################
# Filter1 and its definitions
#######################################################
lcf = open(filterdir+filter1+'.txt','r')      # defintion of the file
riga = lcf.readlines()             # list of lines
riga1 = riga[4:len(riga)]  #list of lines where the wave and transmission are stored
lcf.close()
zp_ef = float(riga[0]) #zero point in energy flux (erg/cm^2/s)
zp_ef_err = zp_ef * 1.0075
filter_ew = riga[1] #equivalent width of the filter
peak_wave = float(riga[2]) # peak wavelength of the filter
system = riga[3] # system used: vega or ab
wavefilter, wavefilter_dmod, transmission= [], [], []
for line in riga1:
    p = line.split()
    wavefilter.append(float(p[0]))
    wavefilter_dmod.append(float(p[0])/(1+float(redshift)))
    transmission.append(float(p[1]))

wavefilterv = array(wavefilter)
wavefilter_dmodv = array(wavefilter_dmod) #plotting purpose
transmissionv = array(transmission)
transmission_initv = array(transmission) #plotting purpose
wavefilter_initv = array(wavefilter) #plotting purpose
fil_obs_min= min(wavefilterv)
fil_obs_max= int(max(wavefilterv)) #integer is needed for a sharper cut-off
#############################################################

band=[1941,2246,2604.57,3561.8,4718.9,6185.2,7499.8,8961.5,3652,4448,5505,6555,7900.4,1524,2320,12370,16471,22126]
filist=['uvw2','uvm2','uvw1','us','gs','rs','is','zs','U','B','V','R','I','FUV','NUV','J','H','K']

print ''
print '##########################################################################'
print 'NB: in a cross K-correction the filter shapes (observed and rest-frame) are the most similar. As a first approximation a wise cross K-correction minimises the difference between the redshifted peak wavelength of the observed filter (value reported in the table below) and the peak wavelength of the rest-frame filter. This programm will evaluate the errors on th K-correction based on those on the redshift, the zero points for the observed and rest-frame filters and the use of a blackbody function instead of a proper spectrum (only if this function is used by the programme). The propagation of the erros is treated as reported in Inserra et al. (2015), ApJ submitted '
print '##########################################################################'
print 'UV GALEX FUV=',band[13], 'NUV=', band[14]
print 'UV Swift+UVOT uvw2=',band[0], 'uvm2=', band[1], 'uvw1=', band[2]
print 'SLOAN  us=',band[3],'gs=',band[4],'rs=',band[5],'is=',band[6],'zs=',band[7]
print 'BESSEL U=',band[8],'B=',band[9],'V=',band[10],'R=',band[11],'I=',band[12]
print 'NIR   J=',band[15],'H=',band[16],'K=',band[17]
print '##########################################################################' 
print ''
print '\033[34mChoose your preferred rest-frame band (cross K-correction is suggested) \033[0m'
print ''

print '\033[34mCentral wavelength rest-frame filter >>> \033[0m', peak_wave/(1+float(redshift))
b = array(band) - peak_wave/(1+redshift)
print '\033[34mFilter suggested for the cross K-correction >>> \033[0m', filist[argmin(abs(b))]
print ''

filter2 = raw_input('To which rest-filter do you want to convert the mangitude FUV,NUV,uvw2,uvm2,uvw1,U,B,V,R,I,us,gs,rs,is,zs,J,H,K [rs] ? ')
if not filter2:
    filter2 = 'rs'

print ''

if filter2 == 'K' or filter2 == 'uvw1' or filter2 == 'uvw2' or filter2 == 'uvm2' or filter2 == 'H' or filter2 =='J':
    questionsys2 = raw_input('Do you want to use the Vega or the AB system ([vega],ab) ? ')
    if not questionsys2:
        questionsys2 = 'vega'

    if questionsys2 == 'ab': #assign the new file with the ZP in AB magnitude
        if filter2 == 'K':
            filter2 = 'K_ab'
        elif filter2 == 'J':
            filter2 = 'J_ab'
        elif filter2 == 'H':
            filter2 = 'H_ab'
        elif filter2 == 'uvw1':
            filter2 = 'uvw1_ab'
        elif filter2 == 'uvw2':
            filter2 = 'uvw2_ab'
        elif filter2 == 'uvm2':
            filter2 ='uvm2_ab'

print ''
#######################################################
# Filter2 and its definitions
#######################################################
lcfr = open(filterdir+filter2+'.txt','r')      # defintion of the file
rigar = lcfr.readlines()             # list of lines
rigar1 = rigar[4:len(rigar)]  #list of lines where wavelength and transmission are stored
lcfr.close()
zp_ef_rest = float(rigar[0]) #zero point in energy flux (erg/cm^2/s)
zp_ef_rest_err = zp_ef_rest * 1.0075
filter_ew_rest = rigar[1] #equivalent width of the filter
peak_wave_rest = float(rigar[2]) # peak wavelength of the filter
system_rest = rigar[3] # system used: vega or ab
wavefilter_rest, transmission_rest= [], []
for line in rigar1:
    p = line.split()
    wavefilter_rest.append(float(p[0]))
    transmission_rest.append(float(p[1]))

wavefilter_restv = array(wavefilter_rest)
transmission_restv = array(transmission_rest)
transmission_restinitv = array(transmission_rest) #plotting purpose
wavefilter_restinitv = array(wavefilter_rest) #plotting purpose
fil_rest_min= min(wavefilter_restv)
fil_rest_max= int(max(wavefilter_restv)) #integer is needed for a sharper cut-off
#############################################################

_ebvg = raw_input('What is the galactic E(B-V) [0.0] ? ')
if not _ebvg:
    ebvg = 0.0
else: 
    ebvg = float(_ebvg)
print ''

_ebvh = raw_input('What is the host E(B-V) [0.0] ? ')
if not _ebvh:
    ebvh = 0.0
else: 
    ebvh = float(_ebvh)
print ''

_sn = raw_input('What SN spectrum ( e.g. 2011ke_) ? ')

#### it recognizes automatically the extension of your file and convert to fits
fileName, fileExtension = os.path.splitext(_sn)
if fileExtension == '.txt' or fileExtension == '.dat' or fileExtension == '.asci' or fileExtension == '.ascii':
    iraf.rspec(_sn,fileName+'.fits',flux='no',dtype='interp')
    sn = fileName+'.fits'
else:
    sn = _sn


#############################
# Magnitude system definitions and check avriable
#############################

print ''
print '#############################################'
print 'Observed frame magnitude system = ',system
print 'Rest frame magnitude system = ',system_rest
print 'Observed passband = ',filter1
print 'Output filter (rest-frame) = ',filter2
print 'SN redshift = ',redshift,'+/-',_redshifterr
print 'E(B-V) galactic = ',ebvg
print 'E(B-V) host = ',ebvh
print '#############################################'
print ''
confort = raw_input('Are you happy with those entries ? ([yes],no) ')
if not confort:
    confort = 'yes'
if confort == 'no' or confort == 'n':
    sys.exit("Bye Bye")
print ''
###################################
##### Mananging the spectrum
###################################

redserrspace1 = redshift+_redshifterr
redserrspace2 = redshift-_redshifterr
spec = sn + "[*,1,1]"            # generally multidimension
iraf.imcopy(sn+'[*,1,1]','sn.fits',verbose='no')             # to create a onedimension fit to use during the script

# redshift correction without absorption                
print '\033[34m*** correcting the spectrum for redshift without Milky Way absorption *** \033[0m'
try:
    iraf.dopcor('sn.fits','sn_dez.fits', redshift=redshift, isveloc='no', flux='no',factor=3)
    iraf.dopcor('sn.fits','sn_dez_err1.fits', redshift=redserrspace1, isveloc='no', flux='no',factor=3)
    iraf.dopcor('sn.fits','sn_dez_err2.fits', redshift=redserrspace2, isveloc='no', flux='no',factor=3)
except:
    print ' WARNING: Problem to redshift the spectrum'

if ebvh == 0.0:
    # galaxy and host reddening correction
    print '\033[31m'+'*** correcting spectrum for galactic reddening ***\033[0m'
    ebv = ebvg+ebvh
    print '\033[31m Total E(B-V) = galactic E(B-V) = \033[0m',ebv
    try:
            iraf.unlearn("deredden")
            iraf.dered('sn_dez.fits',"sn_dez_dered.fits", value=ebv, R=3.1, type='E(B-V)',overrid='yes',uncorre='no')
            iraf.dered('sn_dez_err1.fits',"sn_dez_dered_err1.fits", value=ebv, R=3.1, type='E(B-V)',overrid='yes',uncorre='no')
            iraf.dered('sn_dez_err2.fits',"sn_dez_dered_err2.fits", value=ebv, R=3.1, type='E(B-V)',overrid='yes',uncorre='no')
    except:
            print 'WARNING: it is not possible to correct the spectrum for galactic reddenning '
            try:
                iraf.unlearn("scopy")
                iraf.scopy('sn_dez.fits',"sn_dez_dered.fits", w1='INDEF', w2='INDEF',format='multispec')
                iraf.scopy('sn_dez_err1.fits',"sn_dez_dered_err1.fits", w1='INDEF', w2='INDEF',format='multispec')
                iraf.scopy('sn_dez_err2.fits',"sn_dez_dered_err2.fits", w1='INDEF', w2='INDEF',format='multispec')
            except:
                print 'WARNING: problem to copy the spectrum or with the spectrum fits format'

else:
    #Galaxy reddening correction
    print '\033[31m'+'*** correcting spectrum for galactic reddening ***\033[0m'
    try:
            iraf.dered(sn + '[*,1,1]',"sn_galdered.fits", value=ebvg, R=3.1, type='E(B-V)')
    except:
            print ' WARNING: it is not possible to correct the spectrum for galactic reddenning '
            try:
                iraf.scopy(sn + '[*,1,1]',"sn_galdered.fits", w1='INDEF', w2='INDEF',format='multispec')
            except:
                print ' WARNING: a problem is appeared to copy the spectrum, problems with the spectrum fits format'
    
    print '\033[34m*** correcting the spectrum for redshift ***\033[0m'
    try:
            iraf.dopcor('sn_galdered.fits','sn_galdered_dez.fits', redshift=redshift, isveloc='no', flux='no',factor=3)
            iraf.dopcor('sn_galdered.fits','sn_galdered_dez_err1.fits', redshift=redserrspace1, isveloc='no', flux='no',factor=3)
            iraf.dopcor('sn_galdered.fits','sn_galdered_dez_err2.fits', redshift=redserrspace2, isveloc='no', flux='no',factor=3)
    except:
            print ' WARNING: Problem to redshift the spectrum'
    
    # host reddening correction
    print '\033[31m*** correcting spectrum for host reddening ***\033[0m'
    ebv = ebvg+ebvh
    print '\033[31m NOW total E(B-V) = \033[0m',ebv
    iraf.hedit("sn_galdered_dez.fits", 'DEREDDEN', add='no', addonly='no', delete='yes', verify ='no', show='no')
    iraf.hedit("sn_galdered_dez_err1.fits", 'DEREDDEN', add='no', addonly='no', delete='yes', verify ='no', show='no')
    iraf.hedit("sn_galdered_dez_err2.fits", 'DEREDDEN', add='no', addonly='no', delete='yes', verify ='no', show='no')
    iraf.dered('sn_galdered_dez.fits',"sn_dez_dered.fits", value=ebvh, R=3.1, type='E(B-V)',overrid='yes')
    iraf.dered('sn_galdered_dez_err1.fits',"sn_dez_dered_err1.fits", value=ebvh, R=3.1, type='E(B-V)',overrid='yes')
    iraf.dered('sn_galdered_dez_err2.fits',"sn_dez_dered_err2.fits", value=ebvh, R=3.1, type='E(B-V)',overrid='yes')

#################################################################
#### Preparation for wavelength check
#################################################################


spectrum=iraf.wspec("sn.fits","sn_xbbody.txt", header='no')
lcf = open('sn_xbbody.txt','r')      
riga = lcf.readlines()             
lcf.close()
wave,flux= [],[]
for line in riga:
    p = line.split()
    wave.append(float(p[0]))
    flux.append(float(p[1]))

wavev = array(wave)
fluxv = array(flux)
waveobs_min= min(wavev)
waveobs_max= max(wavev)

spectrum=iraf.wspec("sn_dez.fits","sn_dez_xbbody.txt", header='no')
lcf = open('sn_dez_xbbody.txt','r')      
riga = lcf.readlines()             
lcf.close()
wave,flux= [],[]
for line in riga:
    p = line.split()
    wave.append(float(p[0]))
    flux.append(float(p[1]))

wavedezv = array(wave)
fluxdezv = array(flux)
waverest_min= min(wavedezv)
waverest_max= max(wavedezv)

spectrum=iraf.wspec("sn_dez_err1.fits","sn_dez_xbbody_err1.txt", header='no')
lcf = open('sn_dez_xbbody_err1.txt','r')      
riga = lcf.readlines()             
lcf.close()
wave,flux= [],[]
for line in riga:
    p = line.split()
    wave.append(float(p[0]))
    flux.append(float(p[1]))

wavedezv_err1 = array(wave)
fluxdezv_err1 = array(flux)

spectrum=iraf.wspec("sn_dez_err2.fits","sn_dez_xbbody_err2.txt", header='no')
lcf = open('sn_dez_xbbody_err1.txt','r')      
riga = lcf.readlines()             
lcf.close()
wave,flux= [],[]
for line in riga:
    p = line.split()
    wave.append(float(p[0]))
    flux.append(float(p[1]))

wavedezv_err2 = array(wave)
fluxdezv_err2 = array(flux)

################################
### Define the different cases for K-correction
################################
split = 0

if ((waveobs_min-fil_obs_min) > 50) or ((fil_obs_max-waveobs_max) > 50):
    print ''
    if (waveobs_min-fil_obs_min) > 50:
        print waveobs_min-fil_obs_min,' Angstrom not covered by the observed spectrum in the blue' 
        waveout = waveobs_min-fil_obs_min
    if (fil_obs_max-waveobs_max) > 50:
        print fil_obs_max-waveobs_max,' Angstrom not covered by the observed spectrum in the red'
        waveout = fil_obs_max-waveobs_max
        ############################################
        # Prevent small exceptions for blue bands or the extreme of the NIR
        ############################################
    if filter1 != 'U' or filter1 != 'us' or filter1 != 'K' or filter1 != 'uvw1' or filter1 != 'uvw2' or filter1 != 'uvm2' or filter1 != 'NUV' or filter1 != 'FUV' or filter1 != 'uvw1_ab' or filter1 != 'uvw2_ab' or filter1 != 'uvm2_ab' or filter1 != 'K_ab':

        ###############################
        ### BBody evaluation of the observed spectrum
        ###############################
        BBparams, covar = curve_fit(bbody,wavev,fluxv,p0=(10000,1E-16)) #intial guess
        T= BBparams[0]
        Area = BBparams[1]
        print '\nBlackbody temperature observed spectrum = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))
        outputname = "bbody_sn_fit.dat" #% T
        file = open(outputname,"w")
        file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
        w,f = [],[]
        for wav in range(900,24005):
            file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
            w.append(wav)
            f.append(bbody(wav,T,Area))

        iraf.rspec('bbody_sn_fit.dat','bbody_sn_fit.fits', title='bbodyfit',flux='no',dtype='interp',crval1=900,cdelt1=1)
        iraf.scombine('bbody_sn_fit.fits,sn.fits,sn.fits,sn.fits', 'bsn_combo.fits',combine='median')
        subplot(211)
        plot(wavev,fluxv,'k-',w,f,'r--')
        ylabel('Flux', size=12)
        subplot(212)
        iraf.wspec('bsn_combo.fits','bsn_combo.txt',header='no')
        lcf = open('bsn_combo.txt','r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga[:len(riga)-5]:
            p = line.split()
            wave.append(float(p[0]))
            flux.append(float(p[1]))

        wavev = array(wave)
        fluxv = array(flux)
        wavesp_min= min(wavev)
        wavesp_max= int(max(wavev)) #needed to avoid problems with interp1d
        plot(wave,flux,'k-')
        ylabel('Flux', size=12)
        print '#######################################'
        print '\033[4mDO THAT IS TOO RISKY!!\033[0m I have created for you a \033[1mbbody_sn_fit.fits\033[0m file with the blackbody fit of your observed spectrum. I would suggest to try with a template covering that wavelength region! However, if you like to gamble we can continue using the spectrum on the bottom plot (NB: You have to close the window to continue).'
        print '#######################################'
        show()
        answer = raw_input('Do you want to continue ? ([yes],no) ')
        if not answer:
            answer = 'yes'
        if answer == 'no' or answer =='n':
            sys.exit("Bye Bye")

        elif answer == 'yes':
            print '#######################################'
            print 'OK, but have in mind that the value could be unreasonable. Thus I would suggest you to also use a template and compare the K-correction values between the two methods.'
            print '#######################################'   
            print '\033[31m\033[7m***BBODY version activated***\033[0m'

        lcf = open('bsn_combo.txt','r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_obs_min and float(line.split()[0]) <= fil_obs_max: #to match the spectrum wavelegnths to those of the filter
                wave.append(float(p[0]))
                flux.append(float(p[1]))
            
        wavev = array(wave)
        fluxv = array(flux)
        wavesp_min= min(wavev)
        wavesp_max= int(max(wavev)) #needed to avoid problems with interp1d

        # interpolating the two responses to match the length and sampling coverage
        conf = conv(wavev,fluxv,wavefilterv,transmissionv,wavesp_min,wavesp_max,fil_obs_min,fil_obs_max)

        ##################################
        ### Evaluating the magnitudes
        ##################################
        flux_obs = max(integrate.cumtrapz(conf[0],conf[1])) # using trapezoidal rule to integrate
        flux_obs_err = flux_obs * (1+(waveout-50)*0.0001)

        print 'Apparent magnitude in observed frame ('+filter1+') = ', -2.5*log10(flux_obs/zp_ef)
        iraf.dopcor('bsn_combo.fits','bsn_combo_dez.fits', redshift=redshift, isveloc='no', flux='no',factor=3)
        iraf.dopcor('bsn_combo.fits','bsn_combo_dez_err1.fits', redshift=redserrspace1, isveloc='no', flux='no',factor=3)
        iraf.dopcor('bsn_combo.fits','bsn_combo_dez_err2.fits', redshift=redserrspace2, isveloc='no', flux='no',factor=3)
        iraf.wspec('bsn_combo_dez.fits','bsn_combo_dez.txt',header='no')
        iraf.wspec('bsn_combo_dez_err1.fits','bsn_combo_dez_err1.txt',header='no')
        iraf.wspec('bsn_combo_dez_err2.fits','bsn_combo_dez_err2.txt',header='no')

        lcf = open('bsn_combo_dez.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

        # interpolating the two responses in the rest wavelength to match the length and sampling coverage
        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

        flux_rest = max(integrate.cumtrapz(conf[0],conf[1]))

        ###### doing that again for the errors  ############################################################
        lcf = open('bsn_combo_dez_err1.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

        flux_rest_err1 = max(integrate.cumtrapz(conf[0],conf[1]))

        lcf = open('bsn_combo_dez_err2.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

        flux_rest_err2 = max(integrate.cumtrapz(conf[0],conf[1]))
        ##################################################################################################

        print 'Apparent magnitude in rest frame ('+filter2+') = ', -2.5*log10(flux_rest/zp_ef_rest)
    
        ##################################
        ### Recap of what parameters have been used
        ##################################

        print "\033[32m"+'Observed passband = ',filter1
        print "\033[32m"+'Output filter (rest-frame)= ',filter2
        print "\033[32m"+'SN redshift = ',redshift,'+/-',_redshifterr
        print '\033[32mBlackbody temperature (observed) = %.0f +\- %.0f K' % (T,np.sqrt(covar[0,0]))
        print "\033[32m"+'NONE reddening has been used in this version\033[0m'
    
        ##################################
        ### K-correction evaluation
        ##################################
        kcorrrest_bb_wor=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)) +(2.5*log10(1+redshift))
        
        kcorrrest_bb_wor_err1=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err1/zp_ef_rest)) +(2.5*log10(1+redserrspace1))
        kcorrrest_bb_wor_err2=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err2/zp_ef_rest)) +(2.5*log10(1+redserrspace2))
        kcorrerror_bb_wor=(abs(kcorrrest_bb_wor_err1 - kcorrrest_bb_wor) + abs(kcorrrest_bb_wor_err2 - kcorrrest_bb_wor))/2
        kcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
        kcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
        Kcorrerr_bb = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)) +(2.5*log10(1+redshift)))-(-2.5*log10(flux_obs_err/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)) +(2.5*log10(1+redshift))))
        kcorrerr = sqrt((kcorrerror_bb_wor**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2 + Kcorrerr_bb**2)/4)

        print '#######################################' 
        print '\033[4m'+'K-correction TOTALLY based on an hybrid (if you want to add this value to the observed mag, you have to flip the sign)\033[0m, without reddening and from observed filter '+filter1+' on the hybrid (SN+observed blackbody) observed spectrum to rest-frame filter '+filter2+' on the hybrid (SN+blackbody) = %0.4f' % (kcorrrest_bb_wor),'+/- %0.4f' % (kcorrerr)
        print '#######################################'
        split = 1

    elif filter1 == 'U' or filter1 == 'us' or filter1 == 'uvw1' or filter1 == 'uvw2' or filter1 == 'uvm2' or filter1 == 'NUV' or filter1 == 'FUV' or filter1 == 'uvw1_ab' or filter1 == 'uvw2_ab' or filter1 == 'uvm2_ab':
        print waveobs_min-fil_obs_min,' Angstrom not covered by the observed spectrum in the blue'
        print ''
        print '#######################################'
        print '\033[4mARE YOU SURE ABOUT THAT?\033[0m You are actually using a filter that does not cover your observed spectrum and is too blue (U, us or UV filter). Use a template spectrum, it is better and safer!!'
        print '#######################################'
        sys.exit("Bye Bye")

    elif filter1 == 'K' or filter1 == 'K_ab':
        print fil_obs_max-waveobs_max,' Angstrom not covered by the observed spectrum in the red'
        print ''
        print '#######################################'
        print '\033[4mARE YOU SURE ABOUT THAT?\033[0m You are actually using a filter that does not cover your observed spectrum and is too red (K). Use a template spectrum, it is better and safer!!'
        print '#######################################'
        sys.exit("Bye Bye")

if split == 0:

    if ((waverest_min-fil_rest_min) > 50) or ((fil_rest_max-waverest_max) > 50):
        print ''
        if (waverest_min-fil_rest_min) > 50:
            print waverest_min-fil_rest_min,' Angstrom not covered by the restframe spectrum in the blue'
            waveout = waverest_min-fil_rest_min
        elif (fil_rest_max-waverest_max) > 50:
            print fil_rest_max-waverest_max,' Angstrom not covered by the restframe spectrum in the red'
            waveout = fil_rest_max-waverest_max
    
        print ''
    
        ###################
        ### BBody of the rest frame spectrum
        ###################
    
        BBparams, covar = curve_fit(bbody,wavedezv,fluxdezv,p0=(10000,1E-16)) #initial guess
        T= BBparams[0]
        Area = BBparams[1]
        print '\nBlackbody temperature = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))
        outputname = "bbody_sn_dez_fit.dat" #% T
        file = open(outputname,"w")
        file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
        w,f = [],[]
        for wav in range(900,24005):
           file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
           w.append(wav)
           f.append(bbody(wav,T,Area))
    
        ######### BBody for the error spectra ################################################
        BBparams, covar = curve_fit(bbody,wavedezv_err1,fluxdezv_err1,p0=(10000,1E-16)) #initial guess
        T= BBparams[0]
        Area = BBparams[1]
        outputname = "bbody_sn_dez_fit_err1.dat" #% T
        file = open(outputname,"w")
        file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
        w,f = [],[]
        for wav in range(900,24005):
           file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
           w.append(wav)
           f.append(bbody(wav,T,Area))

        BBparams, covar = curve_fit(bbody,wavedezv_err2,fluxdezv_err2,p0=(10000,1E-16)) #initial guess
        T= BBparams[0]
        Area = BBparams[1]
        outputname = "bbody_sn_dez_fit_err2.dat" #% T
        file = open(outputname,"w")
        file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
        w,f = [],[]
        for wav in range(900,24005):
           file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
           w.append(wav)
           f.append(bbody(wav,T,Area))
        ########################################################################################   

        iraf.rspec('bbody_sn_dez_fit.dat','bbody_sn_dez_fit.fits', title='bbodyfit',flux='no',dtype='interp',crval1=900,cdelt1=1)
        iraf.rspec('bbody_sn_dez_fit_err1.dat','bbody_sn_dez_fit_err1.fits', title='bbodyfit',flux='no',dtype='interp',crval1=900,cdelt1=1)
        iraf.rspec('bbody_sn_dez_fit_err2.dat','bbody_sn_dez_fit_err2.fits', title='bbodyfit',flux='no',dtype='interp',crval1=900,cdelt1=1)
        iraf.scombine('bbody_sn_dez_fit.fits,sn_dez.fits,sn_dez.fits,sn_dez.fits', 'bsn_combo_dez.fits',combine='median')
        iraf.scombine('bbody_sn_dez_fit_err1.fits,sn_dez_err1.fits,sn_dez_err1.fits,sn_dez_err1.fits', 'bsn_combo_dez_err1.fits',combine='median')
        iraf.scombine('bbody_sn_dez_fit_err2.fits,sn_dez_err2.fits,sn_dez_err2.fits,sn_dez_err2.fits', 'bsn_combo_dez_err2.fits',combine='median')
        subplot(211)
        plot(wavedezv,fluxdezv,'k-',w,f,'r--')
        ylabel('Flux', size=12)
        subplot(212)
        iraf.wspec('bsn_combo_dez.fits','bsn_combo_dez.txt',header='no')
        iraf.wspec('bsn_combo_dez_err1.fits','bsn_combo_dez_err1.txt',header='no')
        iraf.wspec('bsn_combo_dez_err2.fits','bsn_combo_dez_err2.txt',header='no')
        lcf = open('bsn_combo_dez.txt','r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga[:len(riga)-5]:
            p = line.split()
            wave.append(float(p[0]))
            flux.append(float(p[1]))
    
        wavev = array(wave)
        fluxv = array(flux)
        plot(wavev,flux,'k-')
        ylabel('Flux', size=12)
        print '#######################################'
        print '\033[4mDO THAT IS TOO RISKY\033[0m (your uncertainties could be bigger than the result)! I have created for you a \033[1mbbody_sn_dez_fit.fits\033[0m file with the blackbody fit of your redshifted (but not dereddened) spectrum. I would suggest to try another iteration with this one! However, if you like to gamble we can continue using the spectrum on the bottom plot (NB: You have to close the window to continue).'
        print '#######################################'
        show()
        answer = raw_input('Do you want to continue ? ([yes],no) ')
        if not answer:
            answer = 'yes'
        if answer == 'no' or answer =='n':
            sys.exit("Bye Bye")
    
        elif answer == 'yes':
            print '#######################################'
            print 'OK, but have in mind that the values could be unreasonable. Thus I would suggest you to also use the blackbody and compare the K-correction values between the two methods.'
            print '#######################################'   
            print '\033[31m\033[7m***BBODY version activated***\033[0m'
    
            ##################################
            ### Evaluating the magnitudes
            ##################################
            iraf.wspec("sn.fits","sn.txt", header='no')
            lcf = open('sn.txt','r')
            riga = lcf.readlines()
            lcf.close()
            wave,flux= [],[]
            for line in riga:
                p = line.split()
                if float(line.split()[0]) >= fil_obs_min and float(line.split()[0]) <= fil_obs_max: #to match the spectrum wavelegnths to those of the filter
                    wave.append(float(p[0]))
                    flux.append(float(p[1]))
                
            wavev = array(wave)
            fluxv = array(flux)
            wavesp_min= min(wavev)
            wavesp_max= int(max(wavev)) #needed to avoid problems with interp1d

            conf = conv(wavev,fluxv,wavefilterv,transmissionv,wavesp_min,wavesp_max,fil_obs_min,fil_obs_max)
    
            flux_obs = max(integrate.cumtrapz(conf[0],conf[1])) # using trapezoidal rule to integrate
            print 'Apparent magnitude in observed frame ('+filter1+') = ', -2.5*log10(flux_obs/zp_ef)

            lcf = open('bsn_combo_dez.txt','r')
            riga = lcf.readlines()
            lcf.close()
            rwave,rflux= [],[]
            for line in riga:
                p = line.split()
                if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                    rwave.append(float(p[0]))
                    rflux.append(float(p[1]))
                
            wave_dezv = array(rwave)
            flux_dezv = array(rflux)
            wavesp_dez_min= min(wave_dezv)
            wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d
    
            # interpolating the two responses in the rest wavelength to match the length and sampling coverage
            conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

            flux_rest = max(integrate.cumtrapz(conf[0],conf[1]))
            flux_rest_err = flux_rest * (1+(waveout-50)*0.0001)
            print 'Apparent magnitude in rest-frame ('+filter2+') = ', -2.5*log10(flux_rest/zp_ef_rest)

            ############ Errors #########################################################################
            lcf = open('bsn_combo_dez_err1.txt','r')
            riga = lcf.readlines()
            lcf.close()
            rwave,rflux= [],[]
            for line in riga:
                p = line.split()
                if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                    rwave.append(float(p[0]))
                    rflux.append(float(p[1]))
                
            wave_dezv = array(rwave)
            flux_dezv = array(rflux)
            wavesp_dez_min= min(wave_dezv)
            wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

            conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

            flux_rest_err1 = max(integrate.cumtrapz(conf[0],conf[1]))

            lcf = open('bsn_combo_dez_err2.txt','r')
            riga = lcf.readlines()
            lcf.close()
            rwave,rflux= [],[]
            for line in riga:
                p = line.split()
                if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                    rwave.append(float(p[0]))
                    rflux.append(float(p[1]))
                
            wave_dezv = array(rwave)
            flux_dezv = array(rflux)
            wavesp_dez_min= min(wave_dezv)
            wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

            conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

            flux_rest_err2 = max(integrate.cumtrapz(conf[0],conf[1]))
            #############################################################################################    
            ##################################
            ### Recap of what parameters have been used
            ##################################
    
            print ''
            print "\033[32m"+'Observed passband = ',filter1
            print "\033[32m"+'Output filter (rest-frame) = ',filter2
            print "\033[32m"+'SN redshift = ',redshift,'+/-',_redshifterr
            print '\033[32mBlackbody temperature = %.0f +\- %.0f K' % (T,np.sqrt(covar[0,0]))
            print "\033[32m"+'NO reddening has been used in this version\033[0m'
            print ''
    
            ##################################
            ### K-correction evaluation
            ##################################
    
            kcorrrest_bb_wor=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))+(2.5*log10(1+redshift))
            kcorrrest_bb_wor_err1=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err1/zp_ef_rest))+(2.5*log10(1+redserrspace1))
            kcorrrest_bb_wor_err2=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err2/zp_ef_rest))+(2.5*log10(1+redserrspace2))
            kcorrerror_bb_wor=(abs(kcorrrest_bb_wor_err1 - kcorrrest_bb_wor) + abs(kcorrrest_bb_wor_err2 - kcorrrest_bb_wor))/2
            kcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
            kcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
            Kcorrerr_bb = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)) +(2.5*log10(1+redshift)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err/zp_ef_rest)) +(2.5*log10(1+redshift))))
            kcorrerr = sqrt((kcorrerror_bb_wor**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2 + Kcorrerr_bb**2)/4)


            print '#######################################' 
            print '\033[4m'+'K-correction based on an hybrid (if you want to add this value to the observed mag, you have to flip the sign)\033[0m, without reddening and from observed filter '+filter1+' on the SN observed spectrum to rest-frame filter '+filter2+' on the hybrid (SN+blackbody) = %0.4f' % (kcorrrest_bb_wor),'+/- %0.4f' % (kcorrerr)
            print '#######################################' 

##################################
### Evaluating the magnitudes
##################################

    else:
        iraf.wspec("sn.fits","sn.txt", header='no')
        lcf = open('sn.txt','r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_obs_min and float(line.split()[0]) <= fil_obs_max: #to match the spectrum wavelegnths to those of the filter
                wave.append(float(p[0]))
                flux.append(float(p[1]))
            
        wavev = array(wave)
        fluxv = array(flux)
        wavesp_min= min(wavev)
        wavesp_max= int(max(wavev)) #needed to avoid problems with interp1d

        conf = conv(wavev,fluxv,wavefilterv,transmissionv,wavesp_min,wavesp_max,fil_obs_min,fil_obs_max)

        flux_obs = max(integrate.cumtrapz(conf[0],conf[1])) # using trapezoidal rule to integrate
        print 'Apparent magnitude in observed frame ('+filter1+') = ', -2.5*log10(flux_obs/zp_ef)
        phot_filtobs_sn = -2.5*log10(flux_obs/zp_ef)

        iraf.wspec("sn_dez_dered.fits","sn_dez_dered.txt", header='no')
        lcf = open('sn_dez_dered.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)

        flux_rest = max(integrate.cumtrapz(conf[0],conf[1]))
        print 'Apparent magnitude in rest-frame ('+filter2+') with reddening applied = ', -2.5*log10(flux_rest/zp_ef_rest)
        phot_filtrest_sn_dez_dered=-2.5*log10(flux_rest/zp_ef_rest)

        #######    Errors  ######################################################################################
        iraf.wspec("sn_dez_dered_err1.fits","sn_dez_dered_err1.txt", header='no')
        lcf = open('sn_dez_dered_err1.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d

        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)
        flux_rest_err1 = max(integrate.cumtrapz(conf[0],conf[1]))

        iraf.wspec("sn_dez_dered_err2.fits","sn_dez_dered_err2.txt", header='no')
        lcf = open('sn_dez_dered_err2.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d
    
        # interpolating the two responses in the rest wavelength to match the length and sampling coverage
        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)
        flux_rest_err2 = max(integrate.cumtrapz(conf[0],conf[1]))
        ######################################################################################################
        phot_filtrest_sn_dez_dered_err1=-2.5*log10(flux_rest_err1/zp_ef_rest)
        phot_filtrest_sn_dez_dered_err2=-2.5*log10(flux_rest_err2/zp_ef_rest)
        
        iraf.wspec("sn_dez.fits","sn_dez.txt", header='no')
        lcf = open('sn_dez.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d
    
        # interpolating the two responses in the rest wavelength to match the length and sampling coverage
        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)
        flux_rest = max(integrate.cumtrapz(conf[0],conf[1]))

        print 'Apparent magnitude in rest-frame ('+filter2+') = ', -2.5*log10(flux_rest/zp_ef_rest)        
        phot_filtrest_sn_dez=-2.5*log10(flux_rest/zp_ef_rest)

        ##################### Errors ######################################################################
        iraf.wspec("sn_dez_err1.fits","sn_dez_err1.txt", header='no')
        lcf = open('sn_dez_err1.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d
    
        # interpolating the two responses in the rest wavelength to match the length and sampling coverage
        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)
        flux_rest_err1 = max(integrate.cumtrapz(conf[0],conf[1]))

        iraf.wspec("sn_dez_err2.fits","sn_dez_err2.txt", header='no')
        lcf = open('sn_dez_err2.txt','r')
        riga = lcf.readlines()
        lcf.close()
        rwave,rflux= [],[]
        for line in riga:
            p = line.split()
            if float(line.split()[0]) >= fil_rest_min and float(line.split()[0]) <= fil_rest_max:
                rwave.append(float(p[0]))
                rflux.append(float(p[1]))
            
        wave_dezv = array(rwave)
        flux_dezv = array(rflux)
        wavesp_dez_min= min(wave_dezv)
        wavesp_dez_max= int(max(wave_dezv)) #needed to avoid problems with interp1d
    
        # interpolating the two responses in the rest wavelength to match the length and sampling coverage
        conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)
        flux_rest_err2 = max(integrate.cumtrapz(conf[0],conf[1]))
        ####################################################################################################
        phot_filtrest_sn_dez_err1=-2.5*log10(flux_rest_err1/zp_ef_rest)
        phot_filtrest_sn_dez_err2=-2.5*log10(flux_rest_err2/zp_ef_rest)
        ##################################
        ### Recap of what parameters have been used
        ##################################
    
        print ''
        print "\033[32m"+'Observed passband = ',filter1
        print "\033[32m"+'Output filter (rest-frame) = ',filter2
        print "\033[32m"+'SN redshift = ',redshift,'+/-',_redshifterr
        print "\033[32m"+'E(B-V) galactic = ',ebvg
        print "\033[32m"+'E(B-V) host = ',ebvh
        print '\033[0m'  
    
        ##################################
        ### Evaluation of K-correction and its errors
        ##################################
    
        kcorrrest_wr=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_dered))+(2.5*log10(1+redshift))
        kcorrrest_wr_err1=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_dered_err1))+(2.5*log10(1+redserrspace1))
        kcorrrest_wr_err2=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_dered_err2))+(2.5*log10(1+redserrspace2))
        kcorrerror_wr=(abs(kcorrrest_wr_err1 - kcorrrest_wr) + abs(kcorrrest_wr_err2 - kcorrrest_wr))/2

        kcorrrest_wor=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez))+(2.5*log10(1+redshift))
        kcorrrest_wor_err1=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_err1))+(2.5*log10(1+redserrspace1))
        kcorrrest_wor_err2=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_err2))+(2.5*log10(1+redserrspace2))
        kcorrerror_wor=(abs(kcorrrest_wor_err1 - kcorrrest_wor) + abs(kcorrrest_wor_err2 - kcorrrest_wor))/2

        kcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
        kcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
        
        kcorrerr_wor = sqrt((kcorrerror_wor**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2)/3)
        kcorrerr_wr = sqrt((kcorrerror_wr**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2)/3)
        
        print '#######################################' 
        print '\033[4m'+'THE REAL K-correction (if you want to add this value to the observed mag, you have to flip the sign) \033[0m is that without reddening from observed filter '+filter1+' to rest-frame filter '+filter2+' = \033[0m %0.4f' % (kcorrrest_wor),'+/- %0.4f' % (kcorrerr_wor)
        print 'K-correction with reddening from observed filter '+filter1+' to rest-frame filter '+filter2+' = %0.4f' % (kcorrrest_wr),'+/- %0.4f' % (kcorrerr_wr)
        print '#######################################'

##################################
### Plotting section
##################################

if (waveobs_min-fil_obs_min) > 50 or (fil_obs_max-waveobs_max) > 50:
     subplot(211)
     listasp= ['bsn_combo.txt','bsn_combo_dez.txt']
     i = 0
     while i != len(listasp):
         wave,flux =[],[]
         lcf = open(listasp[i],'r')
         riga = lcf.readlines()
         lcf.close()
         wave,flux= [],[]
         for line in riga[:len(riga)-5]:
             p = line.split()
             wave.append(float(p[0]))
             flux.append(float(p[1]))
         wavev = array(wave)
         fluxv = array(flux)
         
         # plotting
         spess = ['1','1']   # thikness plot
         tintatipo =['-','-']
         col =['grey','brown']
         plot(wavev,fluxv,color=col[i], ls=tintatipo[i],lw=float(spess[i]))
         legend(['observed','redshifted'],ncol=1,numpoints=1)
         i = i+1
     
     ylabel('Flux', size=14)
     #plot bands and normalized spectra
     subplot(212)
     listasp= ['bsn_combo.txt','bsn_combo_dez.txt']
     i = 0
     while i != len(listasp):
         wave,flux =[],[]
         lcf = open(listasp[i],'r')
         riga = lcf.readlines()
         lcf.close()
         wave,flux= [],[]
         for line in riga[:len(riga)-5]:
             p = line.split()
             wave.append(float(p[0]))
             flux.append(float(p[1]))
         wavev = array(wave)
         fluxv = array(flux)
         normflux = fluxv/max(fluxv)
         # plotting

         spess = ['1','1']   # thikness plot
         tintatipo =['-','-']
         col =['grey','brown']
         plot(wavev,normflux,color=col[i], ls=tintatipo[i],lw=float(spess[i]))
         legend(['observed','redshifted'],ncol=1,numpoints=1)
         i = i+1
    
     #  observed pass-band
     
     plot(wavefilter_initv,transmission_initv,'k:')
     fill_between(wavefilter_initv,transmission_initv,0,where=None,alpha=0.5,color='g')
     
     #  observed pass-band after redshift
     
     plot(wavefilter_dmodv,transmission_initv,'k:')
     fill_between(wavefilter_dmodv,transmission_initv,0,where=None,alpha=0.5,color='yellow')        
     
     #  restframe pass-band

     plot(wavefilter_restinitv,transmission_restinitv,'k:')
     fill_between(wavefilter_restinitv,transmission_restinitv,0,where=None,alpha=0.5,color='c')  

     xlabel('Green ='+filter1+' obs filter, Yellow ='+filter1+' redshifted obs filter, Cyan ='+filter2+' rest filter', size=14)
     ylabel('Normalized Flux', size=14)

elif ((waverest_min-fil_rest_min) > 50) or ((fil_rest_max-waverest_max) > 50):
    subplot(211)
    listasp= ['sn.txt','bsn_combo_dez.txt']
    i = 0
    while i != len(listasp):
        wave,flux =[],[]
        lcf = open(listasp[i],'r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga[:len(riga)-5]:
            p = line.split()
            wave.append(float(p[0]))
            flux.append(float(p[1]))
        wavev = array(wave)
        fluxv = array(flux)
        
        # plotting
        spess = ['1','1']   # thikness plot
        tintatipo =['-','-'] #style
        col =['k','orange'] #colours
        plot(wavev,fluxv,color=col[i], ls=tintatipo[i],lw=float(spess[i])) #plot
        legend(['observed','redshifted'],ncol=1,numpoints=1) #legend
        i = i+1

    ylabel('Flux', size=14)
    #plot bands and normalized spectra
    subplot(212)
    listasp= ['sn.txt','bsn_combo_dez.txt']
    i = 0
    while i != len(listasp):
        wave,flux =[],[]
        lcf = open(listasp[i],'r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga[:len(riga)-5]:
            p = line.split()
            wave.append(float(p[0]))
            flux.append(float(p[1]))
        wavev = array(wave)
        fluxv = array(flux)
        normflux = fluxv/max(fluxv)
        # plotting
        spess = ['1','1']   # thikness plot
        tintatipo =['-','-'] #style
        col =['k','orange'] #colours
        plot(wavev,normflux,color=col[i], ls=tintatipo[i],lw=float(spess[i])) #plot
        legend(['observed','redshifted'],ncol=1,numpoints=1) #legend
        i = i+1

    #  observed pass-band 

    plot(wavefilter_initv,transmission_initv,'k:')
    fill_between(wavefilter_initv,transmission_initv,0,where=None,alpha=0.5,color='g')    

    #  observed pass-band after redshift

    plot(wavefilter_dmodv,transmission_initv,'k:')
    fill_between(wavefilter_dmodv,transmission_initv,0,where=None,alpha=0.5,color='yellow')       
    
    #  restframe pass-band
    
    plot(wavefilter_restinitv,transmission_restinitv,'k:')
    fill_between(wavefilter_restinitv,transmission_restinitv,0,where=None,alpha=0.5,color='c')  
    xlabel('Green ='+filter1+' obs filter, Yellow ='+filter1+' redshifted obs filter, Cyan ='+filter2+' rest filter', size=14)
    ylabel('Normalized Flux', size=14)

    
else:
    subplot(211)
    listasp= ['sn.txt','sn_dez_dered.txt','sn_dez.txt']
    i = 0
    while i != len(listasp):
        wave,flux =[],[]
        lcf = open(listasp[i],'r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga:
            p = line.split()
            wave.append(float(p[0]))
            flux.append(float(p[1]))
        wavev = array(wave)
        fluxv = array(flux)
    
        spess = ['1','1','1']   # thikness plot
        tintatipo =['k-','r-','b-']
        plot(wavev,fluxv, tintatipo[i],lw=float(spess[i]))
        legend(['observed','redshifted+dered','redshifted'],ncol=1,numpoints=1)
        i = i+1

    ylabel('Flux', size=14)
    #plot bands and normalized spectra
    subplot(212)
    listasp= ['sn.txt','sn_dez_dered.txt','sn_dez.txt']
    i = 0
    while i != len(listasp):
        wave,flux =[],[]
        lcf = open(listasp[i],'r')
        riga = lcf.readlines()
        lcf.close()
        wave,flux= [],[]
        for line in riga:
            p = line.split()
            wave.append(float(p[0]))
            flux.append(float(p[1]))
        wavev = array(wave)
        fluxv = array(flux)
        normflux = fluxv/max(fluxv)
        # plotting

        spess = ['1','1','1']   # thikness plot
        tintatipo =['k-','r-','b-']
        plot(wavev,normflux, tintatipo[i],lw=float(spess[i]))
        legend(['observed','redshifted+dered','redshifted'],ncol=1,numpoints=1)
        i = i+1
    
    #  observed pass-band 

    plot(wavefilter_initv,transmission_initv,'k:')
    fill_between(wavefilter_initv,transmission_initv,0,where=None,alpha=0.5,color='g')      
    
    #  observed pass-band after redshift

    plot(wavefilter_dmodv,transmission_initv,'k:')
    fill_between(wavefilter_dmodv,transmission_initv,0,where=None,alpha=0.5,color='yellow')      
    
    #  restframe pass-band

    plot(wavefilter_restinitv,transmission_restinitv,'k:')
    fill_between(wavefilter_restinitv,transmission_restinitv,0,where=None,alpha=0.5,color='c')      
    xlabel('Green ='+filter1+' obs filter, Yellow ='+filter1+' redshifted obs filter, Cyan ='+filter2+' rest filter', size=14)
    ylabel('Normalized Flux', size=14)


#second cleaning process (necessary only for the first time)
os.system('rm -rf sn.txt')
os.system('rm -rf sn.fits')
os.system('rm -rf sn_xbbody.txt')
os.system('rm -rf sn_dez_xbbody.txt')
os.system('rm -rf sn_dez_xbbody_err1.txt')
os.system('rm -rf sn_dez_xbbody_err2.txt')
os.system('rm -rf bbody_sn_dez_fit.dat')
os.system('rm -rf bbody_sn_dez_fit_err1.dat')
os.system('rm -rf bbody_sn_dez_fit_err2.dat')
os.system('rm -rf bbody_sn_dez_fit.fits')
os.system('rm -rf bbody_sn_dez_fit_err1.fits')
os.system('rm -rf bbody_sn_dez_fit_err2.fits')
os.system('rm -rf bbody_sn_fit.fits')
os.system('rm -rf bbody_sn_fit.dat')
os.system('rm -rf sn_dez.txt')
os.system('rm -rf sn_dez_err1.txt')
os.system('rm -rf sn_dez_err2.txt')
os.system('rm -rf sn_dez.fits')
os.system('rm -rf sn_dez_err1.fits')
os.system('rm -rf sn_dez_err2.fits')
os.system('rm -rf sn_dez_dered.fits')
os.system('rm -rf sn_dez_dered_err1.fits')
os.system('rm -rf sn_dez_dered_err2.fits')
os.system('rm -rf sn_dez_dered.txt')
os.system('rm -rf sn_dez_dered_err1.txt')
os.system('rm -rf sn_dez_dered_err2.txt')
os.system('rm -rf sn_galdered_dez.fits')
os.system('rm -rf sn_galdered_dez_err1.fits')
os.system('rm -rf sn_galdered_dez_err2.fits')
os.system('rm -rf sn_galdered.fits')
os.system('rm -rf bsn_combo_dez.fits')
os.system('rm -rf bsn_combo_dez_err1.fits')
os.system('rm -rf bsn_combo_dez_err2.fits')
os.system('rm -rf bsn_combo_dez.txt')
os.system('rm -rf bsn_combo_dez_err1.txt')
os.system('rm -rf bsn_combo_dez_err2.txt')
os.system('rm -rf bsn_combo.fits')
os.system('rm -rf bsn_combo.txt')

show()
