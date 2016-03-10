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
from matplotlib.font_manager import FontProperties
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

description = " P-correction for flux calibrated spectra. List of files are accepted "
usage = "%prog "
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description, version="%prog " + str(s3.__version__))
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    parser.add_option("-s", "--sleep",dest="sleepc", action="store", type="float", default=None, 
                  help='Change the sleep time between cycles. Default is 1s (good for 4GB of RAM or greater), the lower your RAM, the higher it should be.')
    option,args = parser.parse_args()

###### moved here because OptionParser --version conflicts with pyraf version########
#what we need from iraf
from pyraf import iraf

########### option to change the python sleep function between cycles #########
if option.sleepc == None:
    _sleepc = 1
else:
    _sleepc = option.sleepc

################ internal description #############

h="######################################################################\n"+\
  "###########  SuperNova Algorithm for Passband-correction   ###########\n"+\
  "##################             S.N.A.P.           ####################\n"+\
  "##########          C. Inserra  v1.1.0  29/10/2015          ##########\n"+\
  "######################################################################\n"+\
  " P-correction based on the formula m(x) = M(y) + DM + P(y,x)\n"+\
  " BE SURE that the spectra are flux calibrated  \n"+ \
  " If you use this code and find it useful, please give a thought \n"+ \
  " to cite it. \n"+ \
  " The reference is Inserra et al. 2015, ApJ submitted \n"+\
  "######################################################################\n"

print h 

####################################################

#the path where the metatabs dat are
filterdir=s3.__path__[0]+'/metadata/' # To set the directory where are the synphot tabs created

# cleaning process
os.system('rm -rf sn.txt')
os.system('rm -rf sn.fits')
os.system('rm -rf sn_xbbody.txt')
os.system('rm -rf sn_xbbody_std.txt')
os.system('rm -rf bbody_sn_std_fit.dat')
os.system('rm -rf bbody_sn_std_fit.fits')
os.system('rm -rf bbody_sn_fit.fits')
os.system('rm -rf bbody_sn_fit.dat')
os.system('rm -rf sn_std.txt')
os.system('rm -rf sn_std.fits')
os.system('rm -rf bsn_combo_std.fits')
os.system('rm -rf bsn_combo_std.txt')
os.system('rm -rf bsn_combo.fits')
os.system('rm -rf bsn_combo.txt')


#######################################################
# Variable definitions
#######################################################
print ''
print '######################################################################'
print 'Instrumental Passbands and filters IDs accepted by the code:'
print 'Bessell = U B V R I'
print 'Sloan = us gs rs is zs ws(only PS1)'
print '######################################################################'
print ''
print '######################################################################'
print 'Telescope IDs accepted by the code:'
print 'NTT TNG LT NOT LCOGT LSQ PS1 SKYMAPPER OGLE ASIAGO'
print '######################################################################'
print ''

question = raw_input('Do you have a list of spectra ? ([yes]/no) ')
if not question:
    question = 'yes'

if question == 'yes' or question == 'y' or question == 'Y' or question == 'Yes' or question == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             
	lcf.close()
	snlist = []
	for line in riga:
		p = line.split()
		snlist.append(p[0])
else:
	files = raw_input('List the spectra to use (space separated list): ')
	snlist = string.split(files)

print ''
questiontel = raw_input('Do you have a list of telescopes ? ([yes]/no) ')
if not questiontel:
    questiontel = 'yes'

if questiontel == 'yes' or questiontel == 'y' or questiontel == 'Y' or questiontel == 'Yes' or questiontel == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             
	lcf.close()
	tel = []
	for line in riga:
		p = line.split()
		tel.append(p[0])
else:
	tellist = raw_input('List the telescopes you want to use (space separated list) or the single telescope that will be used for all the spectra: ')
	tellist_1 = string.split(tellist)
	if len(tellist_1) != len(snlist):
		if len(tellist_1) == 1:
			tel = tellist_1 * len(snlist)
	else:
		tel = tellist_1
print ''
questionfilobs = raw_input('Do you have a list of observed filters ? ([yes]/no) ')
if not questionfilobs:
    questionfilobs = 'yes'

if questionfilobs == 'yes' or questionfilobs == 'y' or questionfilobs == 'Y' or questionfilobs == 'Yes' or questionfilobs == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()             
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
questionfilrest = raw_input('Do you have a list of standard instrumental passbands ? ([yes]/no) ')
if not questionfilrest:
    questionfilrest = 'yes'

if questionfilrest == 'yes' or questionfilrest == 'y' or questionfilrest == 'Y' or questionfilrest == 'Yes' or questionfilrest == 'YES':
	files = raw_input('List ? [e.g. list, list.txt, list.dat] ')
	lcf = open(files,'r')  
	riga = lcf.readlines()            
	lcf.close()
	frest = []
	for line in riga:
		p = line.split()
		frest.append(p[0])
else:
	frlist = raw_input('List the standard instrumental passbands you want to use (space separated list) or the single passband that will be used for all the spectra: ')
	frlist_1 = string.split(frlist)
	if len(frlist_1) != len(snlist):
		if len(frlist_1) == 1:
			frest = frlist_1 * len(snlist)
	else:
		frest = frlist_1

###########################
### Initialise the arrays
###########################

length = shape(snlist)[0]
pcorr = array(zeros(length))
pcorr_e = array(zeros(length))
method = [None] * len(snlist)
btemp = [None] * len(snlist)
anguncov = [None] * len(snlist)
uncovside = [None] * len(snlist)


##########################
### Creating a txt file
#########################
Tnow = int(strftime("%H%M%S"))
Tnowd = int(strftime("%d%m%Y"))
kcf = "Pcorrection_%.0i_%.0i.txt" % (Tnowd,Tnow)
filekc = open(kcf,"w")
filekc.write("# P-corrections from telescope passband to standard instrumental passbands \n")
filekc.write("# File\tTelescope\tObs-filter\tSTD-passband\tP-corr\terror\tSNAP mode\tBlackbody Temperature\t Angstroms uncovered in the wavelength region\n\n")

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

	filter1 = fobs[ii]
	filter2 = frest[ii]
	tels = tel[ii]
	snum = 1+ii
	print ''
	print '\033[1mSpectrum number\033[0m ', 1+ii
	print 'Spectrum = ', snname, ' & Telescope = ', tels
	sn = snname

	############################# Safety loop to check again if you have everything removed and avoid errors in the programme ###################
	filetoremove = ['sn.txt','sn.fits','sn_xbbody.txt','sn_xbbody_std.txt','bbody_sn_std_fit.dat','bbody_sn_std_fit.fits','bbody_sn_fit.fits','bbody_sn_fit.dat','sn_std.txt','sn_std.fits','bsn_combo_std.fits','bsn_combo_std.txt','bsn_combo.fits','bsn_combo.txt']
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
	lcf = open(filterdir+tels+'/'+filter1+'.txt','r')      # defintion of the file
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

	###################################
	##### Mananging the spectra
	#############################################################

	spec = sn + "[*,1,1]"            # generally multidimension
	iraf.imcopy(sn+'[*,1,1]','sn.fits',verbose='no')             # to create a onedimension fit to use during the script            
	
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

	spectrum=iraf.wspec("sn.fits","sn_xbbody_std.txt", header='no')	
	lcf = open('sn_xbbody_std.txt','r')      
	riga = lcf.readlines()             
	lcf.close()
	wavedez,fluxdez= [],[]
	for line in riga:
		p = line.split()
		wavedez.append(float(p[0]))
		fluxdez.append(float(p[1]))
	
	wavevdez = array(wavedez)
	fluxvdez = array(fluxdez)
	waverest_min= min(wavevdez)
	waverest_max= max(wavevdez)
	
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
		    print '\nBlackbody temperature for telescope filter = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))
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

		    iraf.wspec('bsn_combo.fits','bsn_combo_std.txt',header='no')
		    lcf = open('bsn_combo_std.txt','r')
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
		    ##################################
		    ### K-correction evaluation
		    ##################################
		    pcorrrest_bb_wor=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))
		    pcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
		    pcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
		    pcorrerr_bb = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs_err/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
		    pcorrerr = sqrt((pcorrerrfilt_obs**2 + pcorrerrfilt_rest**2 + pcorrerr_bb**2)/3)

		    pcorr[ii] = pcorrrest_bb_wor
		    pcorr_e[ii] = pcorrerr
		    method[ii] = 'Hybrid Tel_Filter spec_BB'
		    btemp[ii] = bbt	    
		    split = 1

	    elif filter1 == 'U' or filter1 == 'u' or filter1 == 'uvw1' or filter1 == 'uvw2' or filter1 == 'uvm2' or filter1 == 'NUV' or filter1 == 'FUV':
		    print waveobs_min-fil_obs_min,' Angstrom not covered by the observed spectrum in the blue'
		    pcorr[ii] = 0.0
		    pcorr_e[ii] = 0.0
		    method[ii] = 'None'
		    btemp[ii] = 'None' 
		    anguncov[ii] = waveobs_min-fil_obs_min
		    uncovside[ii] = 'Blue Tel-Filter'
	    elif filter1 == 'K':
		    print fil_obs_max-waveobs_max,' Angstrom not covered by the observed spectrum in the red'
		    pcorr[ii] = 0.0
		    pcorr_e[ii] = 0.0
		    method[ii] = 'None'
		    btemp[ii] = 'None'
		    anguncov[ii] = fil_obs_max-waveobs_max
		    uncovside[ii] = 'Red Tel-Filter'

	if split == 0:
		if ((waverest_min-fil_rest_min) > 50) or ((fil_rest_max-waverest_max) > 50):
			print ''
			if (waverest_min-fil_rest_min) > 50:
				print waverest_min-fil_rest_min,' Angstrom not covered by the spectrum in the blue'
				anguncov[ii] = waverest_min-fil_rest_min
				uncovside[ii] = 'Blue STD-Passband'
			elif (fil_rest_max-waverest_max) > 50:
				print fil_rest_max-waverest_max,' Angstrom not covered by the spectrum in the red'
				anguncov[ii] = fil_rest_max-waverest_max
				uncovside[ii] = 'Red STD-Passband'

			print ''
			BBparams, covar = curve_fit(bbody,wavevdez,fluxvdez,p0=(10000,1E-16))
			T= BBparams[0]
			Area = BBparams[1]
			print '\nBlackbody temperature spectrum for standard passband = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))
			bbt = 'BBstd = %.0f +\- %.0f K' % (T,np.sqrt(covar[0,0]))
			outputname = "bbody_sn_std_fit.dat" #% T
			file = open(outputname,"w")
			file.write("# Blackbody temperature = %.0f +\- %.0f K\n" % (T,np.sqrt(covar[0,0])))
			w,f = [],[]
			for wav in range(900,26000):
				file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
				w.append(wav)
				f.append(bbody(wav,T,Area))
		
			iraf.rspec('bbody_sn_std_fit.dat','bbody_sn_std_fit.fits', title='bbodyfit',flux='no',dtype='interp',crval1=900,cdelt1=1)
			iraf.scombine('bbody_sn_std_fit.fits,sn.fits,sn.fits,sn.fits', 'bsn_combo_std.fits',combine='median')
			iraf.wspec('bsn_combo_std.fits','bsn_combo_std.txt',header='no')

		
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

			lcf = open('bsn_combo_std.txt','r')
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
			##################################
			### P-correction evaluation
			##################################
			pcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			pcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			pcorrerr_bb = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest_err/zp_ef_rest))))
			pcorrerr = sqrt((pcorrerrfilt_obs**2 + pcorrerrfilt_rest**2 + pcorrerr_bb**2)/3)
			pcorrrest_bb_wor=-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))
			
			pcorr[ii] = pcorrrest_bb_wor
			pcorr_e[ii] = pcorrerr
			method[ii] = 'Hybrid STD-Passband spec_BB'
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

			iraf.wspec("sn.fits","sn_std.txt", header='no')
			lcf = open('sn_std.txt','r')
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
			phot_filtrest_sn_stdpass=-2.5*log10(flux_rest/zp_ef_rest)
			
			##################################
			### P-correction evaluation
			##################################
			
			pcorrrest_wor=(float(phot_filtobs_sn)-float(phot_filtrest_sn_stdpass))

			pcorrerrfilt_obs = abs((-2.5*log10(flux_obs/zp_ef_err) - (-2.5*log10(flux_rest/zp_ef_rest)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			pcorrerrfilt_rest = abs((-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest_err)))-(-2.5*log10(flux_obs/zp_ef) - (-2.5*log10(flux_rest/zp_ef_rest))))
			pcorrerr_wor = sqrt((pcorrerrfilt_obs**2 + pcorrerrfilt_rest**2)/2)

			pcorr[ii] = pcorrrest_wor
			pcorr_e[ii] = pcorrerr_wor
			method[ii] = 'specTOspec'
			btemp[ii] = 'None'
			anguncov[ii] = 0.0
			uncovside[ii] = 'None'

	# cleaning process (to avoid any problems with small RAM)
	os.system('rm -rf sn.txt')
	os.system('rm -rf sn.fits')
	os.system('rm -rf sn_xbbody.txt')
	os.system('rm -rf sn_xbbody_std.txt')
	os.system('rm -rf bbody_sn_std_fit.dat')
	os.system('rm -rf bbody_sn_std_fit.fits')
	os.system('rm -rf bbody_sn_fit.fits')
	os.system('rm -rf bbody_sn_fit.dat')
	os.system('rm -rf sn_std.txt')
	os.system('rm -rf sn_std.fits')
	os.system('rm -rf bsn_combo_std.fits')
	os.system('rm -rf bsn_combo_std.txt')
	os.system('rm -rf bsn_combo.fits')
	os.system('rm -rf bsn_combo.txt')
	sleep(_sleepc) #to avoid missing files and correction due to a combo of two different spectra (a.k.a. the code arrives at the right step before the system remove the file)
	##########################
	### Adding values to the txt file
	#########################
	filekc.write(snname)
	filekc.write("\t")
	filekc.write(tels)
	filekc.write("\t\t")
	filekc.write(filter1)
	filekc.write("\t\t")
	filekc.write(filter2)
	filekc.write("\t\t")
	filekc.write("%0.3f" % (pcorr[ii]))
	filekc.write("\t\t")
	filekc.write("%0.3f" % (pcorr_e[ii]))
	filekc.write("\t\t")
	filekc.write(method[ii])
	filekc.write("\t\t")
	filekc.write(btemp[ii])
	filekc.write("\t\t")
	filekc.write("%s" % (anguncov[ii]))
	filekc.write("\t")
	filekc.write(uncovside[ii])
	filekc.write("\n")
	###############################
	### plotting section
	###############################
	ax = axes([0.1, 0.1, 0.65, 0.80])
	plot(9999,9999,color='k',marker='o',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='s',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='^',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='v',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='h',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='D',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='d',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='H',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='>',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='<',markeredgecolor='k',ms=10,ls='None')
	snumbl,pcors = [],[]
	if frest[ii] == 'rs':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tels == 'NTT':
			plot(snumbl,pcors,color='orange',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='orange',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='orange',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='orange',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='orange',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='orange',marker='D',ms=12)
		if tel[ii] == 'LSQ':
			plot(snumbl,pcors,color='orange',marker='d',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='orange',marker='H',ms=12)
		if tel[ii] == 'SKYMAPPER':
			plot(snumbl,pcors,color='orange',marker='<',ms=12)

	if frest[ii] == 'is':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tels == 'NTT':
			plot(snumbl,pcors,color='r',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='r',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='r',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='r',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='r',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='r',marker='D',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='r',marker='H',ms=12)
		if tel[ii] == 'SKYMAPPER':
			plot(snumbl,pcors,color='r',marker='<',ms=12)

	if frest[ii] == 'zs':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='#493D26',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='#493D26',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='#493D26',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='#493D26',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='#493D26',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='#493D26',marker='D',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='#493D26',marker='H',ms=12)
		if tel[ii] == 'SKYMAPPER':
			plot(snumbl,pcors,color='#493D26',marker='<',ms=12)

	if frest[ii] == 'gs':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='g',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='g',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='g',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='g',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='g',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='g',marker='D',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='g',marker='H',ms=12)
		if tel[ii] == 'SKYMAPPER':
			plot(snumbl,pcors,color='g',marker='<',ms=12)

	if frest[ii] == 'V':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='#FDD017',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='#FDD017',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='#FDD017',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='#FDD017',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='#FDD017',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='#FDD017',marker='D',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='#FDD017',marker='H',ms=12)
		if tel[ii] == 'OGLE':
			plot(snumbl,pcors,color='#FDD017',marker='>',ms=12)

	if frest[ii] == 'I':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='m',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='m',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='m',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='m',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='m',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='m',marker='D',ms=12)
		if tel[ii] == 'OGLE':
			plot(snumbl,pcors,color='m',marker='>',ms=12)

	if frest[ii] == 'B':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='b',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='b',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='b',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='b',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='b',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='b',marker='D',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='b',marker='H',ms=12)

	if frest[ii] == 'U':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='purple',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='purple',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='purple',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='purple',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='purple',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='purple',marker='D',ms=12)
		if tel[ii] == 'ASIAGO':
			plot(snumbl,pcors,color='purple',marker='H',ms=12)

	if frest[ii] == 'us':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'NTT':
			plot(snumbl,pcors,color='cyan',marker='o',ms=12)
		if tel[ii] == 'LT':
			plot(snumbl,pcors,color='cyan',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='cyan',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='cyan',marker='V',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='cyan',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='cyan',marker='D',ms=12)
		if tel[ii] == 'SKYMAPPER':
			plot(snumbl,pcors,color='cyan',marker='<',ms=12)

	if frest[ii] == 'R':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tels == 'NTT':
			plot(snumbl,pcors,color='#C35817',marker='o',ms=12)
		if tel[ii] == 'LT':
			splot(snumbl,pcors,color='#C35817',marker='s',ms=12)
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='#C35817',marker='^',ms=12)
		if tel[ii] == 'TNG':
			plot(snumbl,pcors,color='#C35817',marker='v',ms=12)
		if tel[ii] == 'LCOGT':
			plot(snumbl,pcors,color='#C35817',marker='h',ms=12)
		if tel[ii] == 'NOT':
			plot(snumbl,pcors,color='#C35817',marker='D',ms=12)

	if frest[ii] == 'ws':
		snumbl.append(snum)
		pcors.append(pcorr[ii])
		if tel[ii] == 'PS1':
			plot(snumbl,pcors,color='grey',marker='^',ms=12)


	ii = ii + 1
then = time.time()
time = then -now
########################
### plotting commands
########################
xl = [0.2,float(len(snlist))+0.1]
yl = [min(pcorr)-0.1,max(pcorr)+0.1]
legend(('NTT', 'LT', 'PS1', 'TNG','LCOGT','NOT','LSQ','ASIAGO','OGLE','SKYMAPPER'), numpoints=1,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
font = FontProperties()
font.set_weight('bold')
text(0.3, max(pcorr)-((max(pcorr)-min(pcorr))/10),'U-Bessell',fontproperties=font,fontsize = 12, color = 'purple')
text(0.3, max(pcorr)-(2*(max(pcorr)-min(pcorr))/10),'B-Bessell',fontproperties=font,fontsize = 12, color = 'b')
text(0.3, max(pcorr)-(3*(max(pcorr)-min(pcorr))/10),'V-Bessell',fontproperties=font,fontsize = 12, color = '#FDD017')
text(0.3, max(pcorr)-(4*(max(pcorr)-min(pcorr))/10),'R-Bessell',fontproperties=font,fontsize = 12, color = '#C35817')
text(0.3, max(pcorr)-(5*(max(pcorr)-min(pcorr))/10),'I-Bessell',fontproperties=font,fontsize = 12, color = 'm')
text(0.3, max(pcorr)-(6*(max(pcorr)-min(pcorr))/10),'u-Sloan',fontproperties=font,fontsize = 12, color = 'c')
text(0.3, max(pcorr)-(7*(max(pcorr)-min(pcorr))/10),'g-Sloan',fontproperties=font,fontsize = 12, color = 'g')
text(0.3, max(pcorr)-(8*(max(pcorr)-min(pcorr))/10),'r-Sloan',fontproperties=font,fontsize = 12, color = 'orange')
text(0.3, max(pcorr)-(9*(max(pcorr)-min(pcorr))/10),'i-Sloan',fontproperties=font,fontsize = 12, color = 'r')
text(0.3, max(pcorr)-(10*(max(pcorr)-min(pcorr))/10),'z-Sloan',fontproperties=font,fontsize = 12, color = '#493D26')
text(0.3, max(pcorr)-(11*(max(pcorr)-min(pcorr))/10),'w-Sloan',fontproperties=font,fontsize = 12, color = 'grey')
xlim(xl[0],xl[1])
ylim(yl[0],yl[1])
xlabel('Spectrum number',size=18)
ylabel('P-correction',size=18)
ax.minorticks_on()
#####################################
##### writing legend on the file
#####################################
filekc.write("\n# -------------------------------------------------------------------------------------------\n")
filekc.write("# Legend for SNAP mode:\n")
filekc.write("# specTOspec                   \t--> P-correction computed with original spectrum for both  \n#\t\t\t\t\t\t            passbands\n")
filekc.write("# Hybrid Tel-Filter spec_BB    \t--> P-correction computed with original spectrum in observed  \n#\t\t\t\t\t\t            frame and SN+Bbody hybrid for the standard passband\n")
filekc.write("# Hybrid STD-Passband spec_BB \t--> P-correction computed with SN+bbody spectrum for observed  \n#\t\t\t\t\t\t            and standard passband\n")
filekc.write("#\n")
filekc.write("# NTT = New Techonology Telescope; LT = Liverpool Telescope; PS1 = Panoramic Survey Telescope and Rapid Response System (Pan-STARSS)  \n# TNG = Telescopio Nazionale Galileo; LCOGT = Las Cumbres Observatory Global Teescope Network\n")
filekc.write("# NOT = Nordic Optical Telescope; LSQ = La SIlla Quest Telescope; ASIAGO = 1.82 Copernico Telescope\n# OGLE = Optical Gravitational Experiment IV; SKYMAPPER = Sky Mapper Telescope \n")
filekc.write("# -------------------------------------------------------------------------------------------\n")

sltime = _sleepc*len(snlist)
print '######################################################################'
print ''
print ' Evaluation done in %.0is, of which %.0is to take a nap to let rest your Random Access Memory (RAM) ' % (time,sltime)	
print ''
print ' \033[46mList of P-corrections\033[0m ' , pcorr
print ' \033[44mVersion used\033[0m ' , method
print ''
print ' A text file has been created ===> Pcorrection_%.0i_%.0i.txt ' % (Tnowd,Tnow)	

show()	

