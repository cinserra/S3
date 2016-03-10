#!/usr/bin/env python

#everything is needed to perform the script and maybe something else
from numpy import *
from scipy import *
from scipy import integrate
from scipy.interpolate import interp1d
import pyfits
import os
import sys
import string 
import shutil
import math
import glob
import s3 #import metadata 
from s3.utilities import *  #import definitions
from time import strftime, sleep
import time
from pylab import *
from scipy.optimize import curve_fit
# pre-set plot parameters, resolution untouched since it is not needed (default=80 dpi) 
from pylab import rcParams
rcParams['figure.figsize'] = 11, 8
rcParams['figure.subplot.top'] = 0.95
rcParams['figure.subplot.right'] = 0.90
rcParams['figure.subplot.left'] = 0.11
###########################################
pypath = os.path.expandvars('$HOME')           # it copies login.cl if it is not in the same dir
if not os.path.isfile('login.cl'):
    shutil.copyfile(pypath+'/iraf/login.cl','login.cl')
###########################################

################### for the help ##################
from optparse import OptionParser

description = " K-correction loop for flux calibrated spectra. List of files are accepted "
usage = "%prog "
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog " + str(s3.__version__))
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    parser.add_option("-s", "--sleep",dest="sleepc", action="store", type="float", default=None, 
                  help='Change the sleep time between cycles. Default is 1s (good for 4GB of RAM or greater), the lower your RAM, the higher it should be.')
    parser.add_option("-r", "--redshifterr",dest="redshifterr", action="store", type="float" ,default=None,
                  help='Change the default error on your redshift (+/- 0.005) to estimate the K-correction errors')
    option,args = parser.parse_args()

###### moved here because OptionParser --version conflicts with pyraf version########
#what we need from iraf
from pyraf import iraf

########### option to change the python sleep function between cycles #########
if option.sleepc == None:
    _sleepc = 1
else:
    _sleepc = option.sleepc

if option.redshifterr == None:
    _redshifterr = 0.005
else:
    _redshifterr = option.redshifterr
################ internal description #############

h="######################################################################\n"+\
  "#########  SuperNova Algorithm for K-correction Evaluation  ##########\n"+\
  "##################     S.N.A.K.E. (loop version)      ################\n"+\
  "##########          C. Inserra  v1.1.0 29/10/2015          ###########\n"+\
  "######################################################################\n"+\
  " K-correction based on the formula m(x) = M(y) + DM + K(y,x)\n"+\
  " BE SURE that the spectra are flux calibrated  \n"+ \
  " If you use this code and find it useful, please give a thought \n"+ \
  " to cite it. \n"+ \
  " The reference is Inserra et al. 2015, ApJ submitted \n"+\
  "######################################################################\n"

print h 

#the path where the metatabs dat are
filterdir=s3.__path__[0]+'/metadata/' # To set the directory

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

question = raw_input('Do you have a list of spectra ? ([yes]/no) ')
if not question:
    question = 'yes'

if question == 'yes' or question == 'y' or question == 'Y' or question == 'Yes' or question == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             # lista di righe intere
	lcf.close()
	snlist = []
	for line in riga:
		p = line.split()
		snlist.append(p[0])
else:
	files = raw_input('List the spectra to use (space separated list): ')
	snlist = string.split(files)

print ''
questionred = raw_input('Do you have a list of redshifts ? ([yes]/no) ')
if not questionred:
    questionred = 'yes'

if questionred == 'yes' or questionred == 'y' or questionred == 'Y' or questionred == 'Yes' or questionred == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             # lista di righe intere
	lcf.close()
	z = []
	for line in riga:
		p = line.split()
		z.append(p[0])
else:
	red = raw_input('List the redshifts you want to use (space separated list) or the single redshift that will be used for all the spectra: ')
	redshift = string.split(red)
	if len(redshift) != len(snlist):
		if len(redshift) == 1:
			z = redshift * len(snlist)
	else:
		z = redshift
		
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
questionfilobs = raw_input('Do you have a list of observed filters ? ([yes]/no) ')
if not questionfilobs:
    questionfilobs = 'yes'

if questionfilobs == 'yes' or questionfilobs == 'y' or questionfilobs == 'Y' or questionfilobs == 'Yes' or questionfilobs == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             # lista di righe intere
	lcf.close()
	fobs = []
	for line in riga:
		p = line.split()
		fobs.append(p[0])
else:
	folist = raw_input('List the observed filters you want to use (space separated list) or the observed filter that will be used for all the spectra: ')
	folist_1 = string.split(folist)
	if len(folist_1) != len(snlist):
		if len(folist_1) == 1:
			fobs = folist_1 * len(snlist)
	else:
		fobs = folist_1
print ''
questionfilrest = raw_input('Do you have a list of rest-frame filters ? ([yes]/no) ')
if not questionfilrest:
    questionfilrest = 'yes'

if questionfilrest == 'yes' or questionfilrest == 'y' or questionfilrest == 'Y' or questionfilrest == 'Yes' or questionfilrest == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             # lista di righe intere
	lcf.close()
	frest = []
	for line in riga:
		p = line.split()
		frest.append(p[0])
else:
	frlist = raw_input('List the rest-frame filters you want to use (space separated list) or the rest-frame filter that will be used for all the spectra: ')
	frlist_1 = string.split(frlist)
	if len(frlist_1) != len(snlist):
		if len(frlist_1) == 1:
			frest = frlist_1 * len(snlist)
	else:
		frest = frlist_1

print ''
print '######################################################################'
print 'Please, bear in mind that the true K-correction should be evaluated   '
print 'without reddening. As a consequence in the cases of K-coorection on   '
print 'hybrid spectra reddening will not be evaluated.'
print '######################################################################'

print ''
_ebvg = raw_input('Galactic E(B-V) [0.0] ? ')
if not _ebvg:
    ebvg = 0.0
else: 
    ebvg = float(_ebvg)
print ''

_ebvh = raw_input('Host E(B-V) [0.0] ? ')
if not _ebvh:
    ebvh = 0.0
else: 
    ebvh = float(_ebvh)
print ''

questionplot = raw_input('Do you want to see your results plotted against redshift ? ([yes]/no) ')
if not questionplot:
    questionplot = 'yes'

###### Arrays definition ##########
length = shape(snlist)[0]
kcorr = array(zeros(length))
kcorr_e = array(zeros(length))
method = [None] * len(snlist)
btemp = [None] * len(snlist)
anguncov = [None] * len(snlist)
uncovside = [None] * len(snlist)
redlist = [None] * len(snlist)


##########################
### Creating a txt file
#########################
Tnow = int(strftime("%H%M%S"))
Tnowd = int(strftime("%d%m%Y"))
kcf = "Kcorrection_%.0i_%.0i.txt" % (Tnowd,Tnow)
filekc = open(kcf,"w")
filekc.write("# K-corrections from observed filter to rest-frame filter \n")
filekc.write("# Galactic E(B-V)=%g, host E(B-V)=%g \n" % (ebvg,ebvh))
filekc.write("# File\tRedshift\tObs-filter\tRest-filter\tK-corr\terror\tSNAKE mode\t Blackbody Temperature\t Angstroms uncovered in the wavelength region\n\n")

now = time.time() 
ii = 0
while ii != len(snlist):
	_snname = snlist[ii]
	#### it recognizes automatically the extension of your file and convert to fits
	fileName, fileExtension = os.path.splitext(_snname)
	if fileExtension == '.txt' or fileExtension == '.dat' or fileExtension == '.asci' or fileExtension == '.ascii':
		iraf.rspec(_snname,fileName+'.fits',flux='no',dtype='interp')
		snname = fileName+'.fits'
	else:
		snname = _snname

	redshift = float(z[ii])
	redserrspace1 = float(redshift)+_redshifterr
	redserrspace2 = float(redshift)-_redshifterr
	filter1 = fobs[ii]
	filter2 = frest[ii]
	print '\033[1mSpectrum number\033[0m ', 1+ii
	print 'Spectrum = ', snname, ' & redshift = ', redshift,'+/-',_redshifterr
 	logterm = 2.5*log10(1+float(redshift))
	sn = snname

	############################# Safety loop to check again if you have everything removed and avoid errors in the programme ###################
	filetoremove = ['sn.txt','sn.fits','sn_xbbody.txt','sn_dez_xbbody.txt','sn_dez_xbbody_err1.txt','sn_dez_xbbody_err2.txt','bbody_sn_dez_fit.dat','bbody_sn_dez_fit_err1.dat','bbody_sn_dez_fit_err2.dat','bbody_sn_dez_fit.fits', \
					'bbody_sn_dez_fit_err1.fits','bbody_sn_dez_fit_err2.fits','bbody_sn_fit.fits','bbody_sn_fit.dat','sn_dez.txt','sn_dez_err1.txt','sn_dez_err2.txt','sn_dez.fits','sn_dez_err1.fits','sn_dez_err2.fits','sn_dez_dered.fits', \
					'sn_dez_dered_err1.fits','sn_dez_dered_err2.fits','sn_dez_dered.txt','sn_dez_dered_err1.txt','sn_dez_dered_err2.txt','sn_galdered_dez.fits','sn_galdered_dez_err1.fits','sn_galdered_dez_err2.fits','sn_galdered.fits', \
					'bsn_combo_dez.fits','bsn_combo_dez_err1.fits','bsn_combo_dez_err2.fits','bsn_combo_dez.txt','bsn_combo_dez_err1.txt','bsn_combo_dez_err2.txt','bsn_combo.fits','bsn_combo.txt']
	jj = 0
	while jj != len(filetoremove):
		if os.path.exists(filetoremove[jj]):
			print ''
			print "######################################################################"
			print "Sorry, I am going too fast for your computer RAM, I need to rest for a bit..."
			print "######################################################################"
			print ''
			for i in xrange(5,0,-1):
				time.sleep(1)
    			sys.stdout.write(str(i)+' ')
    			sys.stdout.flush()
			if os.path.exists(filetoremove[jj]):
				print ''
				print "######################################################################"
				print "Ooops, that is kind of embarassing, apparently there is this file "+filetoremove[jj]+" that is delaying my job. May I ask you to assist me and remove it?"
				for i in xrange(10,0,-1):
					time.sleep(1)
    				sys.stdout.write(str(i)+' ')
    				sys.stdout.flush()

				print "######################################################################"
				print ''
		jj = jj + 1
	#########################################################################################################
	
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
	peak_wave = float(riga[2]) #peak wavelength of the filter
	system = riga[3] # system used: vega or ab
	wavefilter, transmission= [], []
	for line in riga1:
	    p = line.split()
	    wavefilter.append(float(p[0]))
	    transmission.append(float(p[1]))
	
	wavefilterv = array(wavefilter)
	transmissionv = array(transmission)
	fil_obs_min= min(wavefilterv)
	fil_obs_max= int(max(wavefilterv)) #integer is needed for a sharper cut-off
	#############################################################
	
	#######################################################
	# Filter2 and its definitions
	#######################################################
	lcfr = open(filterdir+filter2+'.txt','r')      # defintion of the file
	rigar = lcfr.readlines()             # list of lines
	rigar1 = rigar[4:len(rigar)]  #list of lines where the wave and transmission are stored
	lcfr.close()
	zp_ef_rest = float(rigar[0]) #zero point in energy flux (erg/cm^2/s)
	zp_ef_rest_err = zp_ef_rest * 1.0075
	filter_ew_rest = rigar[1] #equivalent width of the filter
	peak_wave_rest = float(rigar[2]) #peak wavelength of the filter
	system_rest = rigar[3] # system used: vega or ab
	wavefilter_rest, transmission_rest= [], []
	for line in rigar1:
	    p = line.split()
	    wavefilter_rest.append(float(p[0]))
	    transmission_rest.append(float(p[1]))
	
	wavefilter_restv = array(wavefilter_rest)
	transmission_restv = array(transmission_rest)
	fil_rest_min= min(wavefilter_restv)
	fil_rest_max= int(max(wavefilter_restv)) #integer is needed for a sharper cut-off
	#############################################################

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
				print 'WARNING: it is not possible to correct the spectrum for salactic reddenning '
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
			print ' WARNING: it is not possible to correct the spectrum for salactic reddenning '
			try:
				iraf.scopy(sn + '[*,1,1]',"sn_galdered.fits", w1='INDEF', w2='INDEF',format='multispec')
			except:
				print ' WARNING: a problem is appeared to copy the spectrum, problems with the spectrum fits format'
		print '\033[34m*** correcting the spectrum for redshift ***\033[0m'
		try:
			iraf.dopcor('sn_galdered.fits','sn_galdered_dez.fits', redshift=redshift, isveloc='no', flux='no')
			iraf.dopcor('sn_galdered.fits','sn_galdered_dez_err1.fits', redshift=redserrspace1, isveloc='no', flux='no',factor=3)
			iraf.dopcor('sn_galdered.fits','sn_galdered_dez_err2.fits', redshift=redserrspace2, isveloc='no', flux='no',factor=3)
		except:
			print ' WARNING: Problem to redshift the spectrum'
		# host reddening correction
		print '\033[31m*** correcting spectrum for host reddening ***\033[0m'
		ebv = ebvg+ebvh
		print '\033[31m NOW total E(B-V) = ',ebv
		iraf.hedit("sn_galdered_dez.fits", 'DEREDDEN', add='no', addonly='no', delete='yes', verify ='no')
		iraf.hedit("sn_galdered_dez_err1.fits", 'DEREDDEN', add='no', addonly='no', delete='yes', verify ='no', show='no')
		iraf.hedit("sn_galdered_dez_err2.fits", 'DEREDDEN', add='no', addonly='no', delete='yes', verify ='no', show='no')
		iraf.dered('sn_galdered_dez.fits',"sn_dez_dered.fits", value=ebvh, R=3.1, type='E(B-V)',overrid='yes')
		iraf.dered('sn_galdered_dez_err1.fits',"sn_dez_dered_err1.fits", value=ebvh, R=3.1, type='E(B-V)',overrid='yes')
		iraf.dered('sn_galdered_dez_err2.fits',"sn_dez_dered_err2.fits", value=ebvh, R=3.1, type='E(B-V)',overrid='yes')
		print '\033[0m' 
	

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
	
	spectrum=iraf.wspec("sn_dez_dered.fits","sn_dez_xbbody.txt", header='no')
	lcf = open('sn_dez_xbbody.txt','r')     
	riga = lcf.readlines()          
	lcf.close()
	wavedez,fluxdez= [],[]
	for line in riga:
		p = line.split()
		wavedez.append(float(p[0]))
		fluxdez.append(float(p[1]))
	
	wavedezv = array(wavedez)
	fluxdezv = array(fluxdez)
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
	split = 0 # splitting value

	if ((waveobs_min-fil_obs_min) > 50) or ((fil_obs_max-waveobs_max) > 50):
	    print ''
	    if (waveobs_min-fil_obs_min) > 50:
	        print waveobs_min-fil_obs_min,' Angstrom not covered by the observed spectrum in the blue' 
	        anguncov[ii] = waveobs_min-fil_obs_min
	        uncovside[ii] = 'Blue'
	    if (fil_obs_max-waveobs_max) > 50:
	        print fil_obs_max-waveobs_max,' Angstrom not covered by the observed spectrum in the red'
	        anguncov[ii] = fil_obs_max-waveobs_max
	        uncovside[ii] = 'Red'
	        ############################################
	        # Prevent small exceptions for blue bands or the extreme of the NIR
	        ############################################
	    if filter1 != 'U' or filter1 != 'u' or filter1 != 'K' or filter1 != 'uvw1' or filter1 != 'uvw2' or filter1 != 'uvm2' or filter1 != 'NUV' or filter1 != 'FUV':

		    ###############################
		    ### BBody evaluation of the observed spectrum
		    ###############################
		    BBparams, covar = curve_fit(bbody,wavev,fluxv,p0=(10000,1E-16)) #intial guess
		    T= BBparams[0]
		    Area = BBparams[1]
		    print '\nBlackbody temperature observed spectrum = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))
		    bbt = 'BBobs = %.0f +\- %.0f K' % (T,np.sqrt(covar[0,0]))
		    outputname = "bbody_sn_fit.dat" #% T
		    file = open(outputname,"w")
		    file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
		    w,f = [],[]
		    for wav in range(900,26000):
		        file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
		        w.append(wav)
		        f.append(bbody(wav,T,Area))

		    iraf.rspec('bbody_sn_fit.dat','bbody_sn_fit.fits', title='bbodyfit',flux='no',dtype='interp',crval1=900,cdelt1=1)
		    iraf.scombine('bbody_sn_fit.fits,sn.fits,sn.fits,sn.fits', 'bsn_combo.fits',combine='median')

		    iraf.wspec('bsn_combo.fits','bsn_combo.txt',header='no')

		    print '#######################################'
		    print '\033[4mSince now you are working with an hybrid spectrum+blackbody, if you want additional information, please use the normal version.\033[0m '
		    print '#######################################'

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
		    flux_obs_err = flux_obs * (1+(anguncov[ii]-50)*0.0001)

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
		    kcorrrest_bb_wor=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)) +logterm
		    kcorrrest_bb_wor_err1=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err1/zp_ef_rest)) +(2.5*log10(1+redserrspace1))
		    kcorrrest_bb_wor_err2=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err2/zp_ef_rest)) +(2.5*log10(1+redserrspace2))
		    kcorrerror_bb_wor=(abs(kcorrrest_bb_wor_err1 - kcorrrest_bb_wor) + abs(kcorrrest_bb_wor_err2 - kcorrrest_bb_wor))/2
		    kcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
		    kcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
		    Kcorrerr_bb = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs_err/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
		    kcorrerr = sqrt((kcorrerror_bb_wor**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2 + Kcorrerr_bb**2)/4)
		    
		    kcorr[ii] = kcorrrest_bb_wor
		    kcorr_e[ii] = kcorrerr
		    method[ii] = 'Hybrid obs_spec_BB'
		    btemp[ii] = bbt		    
		    split = 1

	    elif filter1 == 'U' or filter1 == 'u' or filter1 == 'uvw1' or filter1 == 'uvw2' or filter1 == 'uvm2' or filter1 == 'NUV' or filter1 == 'FUV':
		    print waveobs_min-fil_obs_min,' Angstrom not covered by the observed spectrum in the blue'
		    kcorr[ii] = 0.0
		    kcorr_e[ii] = 0.0
		    method[ii] = 'None'
		    btemp[ii] = 'None' 
		    anguncov[ii] = waveobs_min-fil_obs_min
		    uncovside[ii] = 'Blue'
	    elif filter1 == 'K':
		    print fil_obs_max-waveobs_max,' Angstrom not covered by the observed spectrum in the red'
		    kcorr[ii] = 0.0
		    kcorr_e[ii] = 0.0
		    method[ii] = 'None'
		    btemp[ii] = 'None'
		    anguncov[ii] = fil_obs_max-waveobs_max
		    uncovside[ii] = 'Red'

	if split == 0:
		if ((waverest_min-fil_rest_min) > 50) or ((fil_rest_max-waverest_max) > 50):
			print ''
			if (waverest_min-fil_rest_min) > 50:
				print waverest_min-fil_rest_min,' Angstrom not covered by the restframe spectrum in the blue'
				anguncov[ii] = waverest_min-fil_rest_min
				uncovside[ii] = 'Blue'
			elif (fil_rest_max-waverest_max) > 50:
				print fil_rest_max-waverest_max,' Angstrom not covered by the restframe spectrum in the red'
				anguncov[ii] = fil_rest_max-waverest_max
				uncovside[ii] = 'Red'

			print ''
			BBparams, covar = curve_fit(bbody,wavedezv,fluxdezv,p0=(10000,1E-16))
			T= BBparams[0]
			Area = BBparams[1]
			print '\nBlackbody temperature rest spectrum = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))
			bbt = 'BBrest = %.0f +\- %.0f K' % (T,np.sqrt(covar[0,0]))
			outputname = "bbody_sn_dez_fit.dat" #% T
			file = open(outputname,"w")
			file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
			w,f = [],[]
			for wav in range(900,26000):
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
			iraf.wspec('bsn_combo_dez.fits','bsn_combo_dez.txt',header='no')
			iraf.wspec('bsn_combo_dez_err1.fits','bsn_combo_dez_err1.txt',header='no')
			iraf.wspec('bsn_combo_dez_err2.fits','bsn_combo_dez_err2.txt',header='no')

			print '#######################################'
			print ' Since now you are using a rest-frame spectrum based on a combination of SN+Blackbody '
			print '#######################################'

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
			
			# interpolating the two responses to match the length and sampling coverage
			conf = conv(wavev,fluxv,wavefilterv,transmissionv,wavesp_min,wavesp_max,fil_obs_min,fil_obs_max)
			flux_obs = max(integrate.cumtrapz(conf[0],conf[1])) # using trapezoidal rule to integrate

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
			flux_rest_err = flux_rest * (1+(anguncov[ii]-50)*0.0001)

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

			kcorrrest_bb_wor=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))+logterm
			kcorrrest_bb_wor_err1=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err1/zp_ef_rest))+(2.5*log10(1+redserrspace1))
			kcorrrest_bb_wor_err2=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err2/zp_ef_rest))+(2.5*log10(1+redserrspace2))
			kcorrerror_bb_wor=(abs(kcorrrest_bb_wor_err1 - kcorrrest_bb_wor) + abs(kcorrrest_bb_wor_err2 - kcorrrest_bb_wor))/2
			kcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			kcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			Kcorrerr_bb = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)) +(2.5*log10(1+redshift)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err/zp_ef_rest)) +(2.5*log10(1+redshift))))
			kcorrerr = sqrt((kcorrerror_bb_wor**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2 + Kcorrerr_bb**2)/4)

			kcorr[ii] = kcorrrest_bb_wor
			kcorr_e[ii] = kcorrerr
			method[ii] = 'Hybrid rest_spec_BB'
			btemp[ii] = bbt
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
			# interpolating the two responses to match the length and sampling coverage

			conf = conv(wavev,fluxv,wavefilterv,transmissionv,wavesp_min,wavesp_max,fil_obs_min,fil_obs_max)

			##################################
			### Evaluating the magnitudes
			##################################
			flux_obs = max(integrate.cumtrapz(conf[0],conf[1])) # using trapezoidal rule to integrate
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
			
			# interpolating the two responses in the rest wavelength to match the length and sampling coverage
			conf = conv(wave_dezv,flux_dezv,wavefilter_restv,transmission_restv,wavesp_dez_min,wavesp_dez_max,fil_rest_min,fil_rest_max)
			flux_rest = max(integrate.cumtrapz(conf[0],conf[1]))

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
 			
 			# interpolating the two responses in the rest wavelength to match the length and sampling coverage
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
			### K-correction evaluation
			##################################
			
			kcorrrest_wr=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_dered))+logterm
			kcorrrest_wr_err1=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_dered_err1))+(2.5*log10(1+redserrspace1))
			kcorrrest_wr_err2=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_dered_err2))+(2.5*log10(1+redserrspace2))
			kcorrerror_wr=(abs(kcorrrest_wr_err1 - kcorrrest_wr) + abs(kcorrrest_wr_err2 - kcorrrest_wr))/2

			kcorrrest_wor=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez))+logterm
			kcorrrest_wor_err1=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_err1))+(2.5*log10(1+redserrspace1))
			kcorrrest_wor_err2=(float(phot_filtobs_sn)-float(phot_filtrest_sn_dez_err2))+(2.5*log10(1+redserrspace2))
			kcorrerror_wor=(abs(kcorrrest_wor_err1 - kcorrrest_wor) + abs(kcorrrest_wor_err2 - kcorrrest_wor))/2
			
			kcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			kcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			
			kcorrerr_wor = sqrt((kcorrerror_wor**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2)/3)
			kcorrerr_wr = sqrt((kcorrerror_wr**2 + kcorrerrfilt_obs**2 + kcorrerrfilt_rest**2)/3)

			kcorr[ii] = kcorrrest_wr
			kcorr_e[ii] = kcorrerr_wr
			method[ii] = 'specTOspec'
			btemp[ii] = 'None'
			anguncov[ii] = 0.0
			uncovside[ii] = 'None'


	# cleaning process (unnecessary)
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
	sleep(_sleepc) #to avoid missing files and correction due to a combo of two different spectra (a.k.a. the code arrives at the right step before the system remove the file)
	##########################
	### Adding values to the txt file
	#########################
	redlist[ii] = redshift
	filekc.write(snname)
	filekc.write("\t")
	filekc.write("%s+/-%s" % (redshift,_redshifterr))
	filekc.write("\t\t")
	filekc.write(filter1)
	filekc.write("\t\t")
	filekc.write(filter2)
	filekc.write("\t\t")
	filekc.write("%0.3f" % (kcorr[ii]))
	filekc.write("\t\t")
	filekc.write("+/-%0.3f" % (kcorr_e[ii]))
	filekc.write("\t\t")
	filekc.write(method[ii])
	filekc.write("\t\t")
	filekc.write(btemp[ii])
	filekc.write("\t\t")
	filekc.write("%s" % (anguncov[ii]))
	filekc.write("\t")
	filekc.write(uncovside[ii])
	filekc.write("\n")
	#ax = subplot(111)
	if questionplot == 'yes':
		ax = axes([0.1, 0.1, 0.65, 0.82])
		plot(9999,9999,color='b',marker='o',markeredgecolor='k',ls='None')
		plot(9999,9999,color='g',marker='d',markeredgecolor='k',ls='None')
		plot(9999,9999,color='r',marker='s',markeredgecolor='k',ls='None')
		reds,redo,redr,kcors,kcoro,kcor1 = [],[],[],[],[],[]
		if method[ii] == 'specTOspec':
			reds.append(redshift)
			kcors.append(kcorr[ii])
			plot(reds,kcors,'bo',ms=12)
		if method[ii] == 'Hybrid rest_spec_BB':		
			redr.append(redshift)
			kcor1.append(kcorr[ii])
			plot(redr,kcor1,'gd',ms=12)
		if method[ii] == 'Hybrid obs_spec_BB':
			redo.append(redshift)
			kcoro.append(kcorr[ii])
			plot(redo,kcoro,'rs',ms=12)
	ii = ii + 1

then = time.time()
time = then -now
################
### plotting commands
##################
if questionplot == 'yes':
	xl = [float(min(redlist))-0.07,float(max(redlist))+0.07]
	yl = [min(kcorr)-0.3,max(kcorr)+0.3]
	legend(('specTOspec', 'Hybrid rest specBB', 'Hybrid obs specBB'), numpoints=1,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
	ltext = gca().get_legend().get_texts()
	setp(ltext[0], fontsize = 14, color = 'b')
	setp(ltext[1], fontsize = 14, color = 'g')
	setp(ltext[2], fontsize = 14, color = 'r')
	xlim(xl[0],xl[1])
	ylim(yl[0],yl[1])
	title('From various observed filters to various rest filters')
	xlabel('Redshift',size=18)
	ylabel('K-correction',size=18)
	ax.minorticks_on()
	show()
####################
##### writing legend on the file
####################
filekc.write("\n# ----------------------------------------------------------------------------------\n")
filekc.write("# Legend for SNAKE mode:\n")
filekc.write("# specTOspec        \t--> K-correction computed with original spectrum in observed  \n#\t\t\t\t\t\tand rest frame\n")
filekc.write("# Hybrid rest_spec_BB\t--> K-correction computed with original spectrum in observed  \n#\t\t\t\t\t\tframe and SN+Bbody hybrid in rest frame\n")
filekc.write("# Hybrid obs_spec_BB \t--> K-correction computed with SN+bbody spectrum in observed  \n#\t\t\t\t\t\tand rest frame\n")
filekc.write("# ----------------------------------------------------------------------------------\n")

sltime = _sleepc*len(snlist)
print '######################################################################'
print ''
print ' Evaluation done in %.0is, of which %.0is to take a nap to let rest your Random Access Memory (RAM) ' % (time,sltime)	
print ''
print ' \033[46mList of K-corrections\033[0m ' , kcorr
print ' \033[44mVersion used\033[0m ' , method
print ''
print ' A text file has been created ===> Kcorrection_%.0i_%.0i.txt ' % (Tnowd,Tnow)	


