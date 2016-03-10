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
from time import strftime, sleep
import time
from pylab import *
from scipy.optimize import curve_fit
import s3 #import metadata 
from s3.utilities import *  #import definitions
# pre-set plot parameters, resolution untouched since it is not needed (default=80 dpi) 
from matplotlib.font_manager import FontProperties
from pylab import rcParams
rcParams['figure.figsize'] = 11, 8
rcParams['figure.subplot.top'] = 0.95
rcParams['figure.subplot.right'] = 0.90
rcParams['figure.subplot.left'] = 0.11
###################################################
pypath = os.path.expandvars('$HOME')           # it copies login.cl if it is not in the same dir
if not os.path.isfile('login.cl'):
    shutil.copyfile(pypath+'/iraf/login.cl','login.cl')
###################################################

################### for the help ##################
from optparse import OptionParser

description = " Synthetic magnitudes from flux calibrated spectra "
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
  "###############  Synthetic  Magnitudes from Spectra  #################\n"+\
  "###################             S.M.S.         #######################\n"+\
  "##########         C. Inserra  v1.1.0  29/10/2015         ############\n"+\
  "######################################################################\n"+\
  " PLEASE READ CAREFULLY                            \n"+ \
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
os.system('rm -rf bbody_sn_fit.fits')
os.system('rm -rf bbody_sn_fit.dat')
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
questionfilobs = raw_input('Do you have a list of filters ? ([yes]/no) ')
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
	folist = raw_input('List the filters you want to use (space separated list) or the observed filter that will be used for all the spectra: ')
	folist_1 = string.split(folist)
	if len(folist_1) != len(snlist):
		if len(folist_1) == 1:
			fobs = folist_1 * len(snlist)
	else:
		fobs = folist_1

length = shape(snlist)[0]
mag = array(zeros(length))
mag_e = array(zeros(length))
method = [None] * len(snlist)
btemp = [None] * len(snlist)
anguncov = [None] * len(snlist)
uncovside = [None] * len(snlist)


##########################
### Creating a txt file
#########################
Tnow = int(strftime("%H%M%S"))
Tnowd = int(strftime("%d%m%Y"))
kcf = "Magnitudes_%.0i_%.0i.txt" % (Tnowd,Tnow)
filekc = open(kcf,"w")
filekc.write("# Synthetic magnitudesfrom spectra \n")
filekc.write("# File\tFilter\tMagnitude\t errore\t SMS mode\t Blackbody Temperature\t Angstroms uncovered in the wavelength region\n\n")

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
	snum = 1+ii
	print ''
	print '\033[1mSpectrum number\033[0m ', 1+ii
	print 'Spectrum = ', snname
	sn = snname

	############################# Safety loop to check again if you have everything removed and avoid errors in the programme ###################
	filetoremove = ['sn.txt','sn.fits','sn_xbbody.txt','bsn_combo.fits','bsn_combo.txt','bbody_sn_fit.dat','bbody_sn_fit.fits']
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
	

	#############################################################
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
		
		    phot_filtobs_bb=-2.5*log10(flux_obs/zp_ef)
		    mcorrerrfilt_obs = abs(-2.5*log10(flux_obs/zp_ef_err) -(-2.5*log10(flux_obs/zp_ef)))
		    mcorrerr_bb = abs(-2.5*log10(flux_obs_err/zp_ef) - (-2.5*log10(flux_obs/zp_ef)))
		    mcorrerr = sqrt((mcorrerrfilt_obs**2 ++ mcorrerr_bb**2)/2)

		    mag[ii] = phot_filtobs_bb
		    mag_e[ii] = mcorrerr
		    method[ii] = 'Hybrid spec_BB'
		    btemp[ii] = bbt	    
		    split = 1

	    elif filter1 == 'U' or filter1 == 'u' or filter1 == 'uvw1' or filter1 == 'uvw2' or filter1 == 'uvm2' or filter1 == 'NUV' or filter1 == 'FUV':
		    print waveobs_min-fil_obs_min,' Angstrom not covered by the observed spectrum in the blue'
		    mag[ii] = 0.0
		    mag_e[ii] = 0.0
		    method[ii] = 'None'
		    btemp[ii] = 'None' 
		    anguncov[ii] = waveobs_min-fil_obs_min
		    uncovside[ii] = 'Blue'
	    elif filter1 == 'K':
		    print fil_obs_max-waveobs_max,' Angstrom not covered by the observed spectrum in the red'
		    mag[ii] = 0.0
		    mag_e[ii] = 0.0
		    method[ii] = 'None'
		    btemp[ii] = 'None'
		    anguncov[ii] = fil_obs_max-waveobs_max
		    uncovside[ii] = 'Red'

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

		mcorrerrfilt_obs = abs(-2.5*log10(flux_obs/zp_ef_err) -(-2.5*log10(flux_obs/zp_ef)))

		mag[ii] = phot_filtobs_sn
		mag_e[ii] = mcorrerrfilt_obs
		method[ii] = 'specTOspec'
		btemp[ii] = 'None'
		anguncov[ii] = 0.0
		uncovside[ii] = 'None'

	# cleaning process (to avoid any problems with small RAM)
	os.system('rm -rf sn.txt')
	os.system('rm -rf sn.fits')
	os.system('rm -rf sn_xbbody.txt')
	os.system('rm -rf bsn_combo.fits')
	os.system('rm -rf bsn_combo.txt')
	os.system('rm -rf bbody_sn_fit.dat')
	os.system('rm -rf bbody_sn_fit.fits')
	sleep(_sleepc) #to avoid missing files and correction due to a combo of two different spectra (a.k.a. the code arrives at the right step before the system remove the file)
	##########################
	### Adding values to the txt file
	#########################
	filekc.write(snname)
	filekc.write("\t")
	filekc.write(filter1)
	filekc.write("\t\t")
	filekc.write("%0.3f" % (mag[ii]))
	filekc.write("\t\t")
	filekc.write("%0.3f" % (mag_e[ii]))
	filekc.write("\t\t")
	filekc.write(method[ii])
	filekc.write("\t\t")
	filekc.write(btemp[ii])
	filekc.write("\t\t")
	filekc.write("%s" % (anguncov[ii]))
	filekc.write("\t")
	filekc.write(uncovside[ii])
	filekc.write("\n")

	ax = axes([0.1, 0.1, 0.65, 0.80])
	plot(9999,9999,color='k',marker='o',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='s',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='^',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='d',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='v',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='D',markeredgecolor='k',ms=10,ls='None')
	plot(9999,9999,color='k',marker='h',markeredgecolor='k',ms=10,ls='None')
	mags,snums = [],[]
	if fobs[ii] == 'rs':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='s',color='orange',ms=12)
	if fobs[ii] == 'is':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='s',color='r',ms=12)
	if fobs[ii] == 'zs':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='s',color='brown',ms=12)
	if fobs[ii] == 'gs':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='s',color='green',ms=12)
	if fobs[ii] == 'us':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='s',color='blue',ms=12)
	if fobs[ii] == 'U':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='o',color='darkblue',ms=12)
	if fobs[ii] == 'B':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='o',color='cyan',ms=12)
	if fobs[ii] == 'V':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='o',color='yellow',ms=12)
	if fobs[ii] == 'R':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='o',color='#C35817',ms=12)
	if fobs[ii] == 'I':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='o',color='m',ms=12)	
	if fobs[ii] == 'J':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='^',color='#6F4E37',ms=12)	
	if fobs[ii] == 'J_ab':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='d',color='#6F4E37',ms=12)	
	if fobs[ii] == 'H':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='^',color='#B87333',ms=12)	
	if fobs[ii] == 'H_ab':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='d',color='#B87333',ms=12)	
	if fobs[ii] == 'K':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='^',color='#827B60',ms=12)	
	if fobs[ii] == 'K_ab':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='d',color='#827B60',ms=12)	
	if fobs[ii] == 'uvw1':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='v',color='#7FFFD4',ms=12)	
	if fobs[ii] == 'uvw1_ab':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='D',color='#7FFFD4',ms=12)	
	if fobs[ii] == 'uvm2':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='v',color='#6960EC',ms=12)	
	if fobs[ii] == 'uvm2_ab':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='D',color='#6960EC',ms=12)	
	if fobs[ii] == 'uvw2':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='v',color='#7D0552',ms=12)	
	if fobs[ii] == 'uvw2_ab':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='D',color='#7D0552',ms=12)	
	if fobs[ii] == 'FUV':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='h',color='#4B0082',ms=12)	
	if fobs[ii] == 'NUV':
		mags.append(mag[ii])
		snums.append(snum)
		plot(snums,mags,marker='h',color='#95B9C7',ms=12)	
	ii = ii + 1

then = time.time()
time = then -now
###########################
### plotting commands
###########################
xl = [0.2,float(len(snlist))+0.1]
yl = [min(mag)-0.6,max(mag)+0.6]
legend(('Bessell', 'Sloan', 'NIR(Vega)', 'NIR(ab)','SwiftUV(Vega)','SwiftUV(ab)','GALEX'), numpoints=1,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
font = FontProperties()
font.set_weight('bold')
text(0.3,max(mag)-((max(mag)-min(mag))/18),'U-Bessell',fontproperties=font,fontsize = 12, color = 'darkblue')
text(0.3,max(mag)-(2*(max(mag)-min(mag))/18),'B-Bessell',fontproperties=font,fontsize = 12, color = 'c')
text(0.3, max(mag)-(3*(max(mag)-min(mag))/18),'V-Bessell',fontproperties=font,fontsize = 12, color = 'yellow')
text(0.3, max(mag)-(4*(max(mag)-min(mag))/18),'R-Bessell',fontproperties=font,fontsize = 12, color = '#C35817')
text(0.3, max(mag)-(5*(max(mag)-min(mag))/18),'I-Bessell',fontproperties=font,fontsize = 12, color = 'm')
text(0.3, max(mag)-(6*(max(mag)-min(mag))/18),'u-Sloan',fontproperties=font,fontsize = 12, color = 'b')
text(0.3, max(mag)-(7*(max(mag)-min(mag))/18),'g-Sloan',fontproperties=font,fontsize = 12, color = 'g')
text(0.3, max(mag)-(8*(max(mag)-min(mag))/18),'r-Sloan',fontproperties=font,fontsize = 12, color = 'orange')
text(0.3, max(mag)-(9*(max(mag)-min(mag))/18),'i-Slaon',fontproperties=font,fontsize = 12, color = 'r')
text(0.3, max(mag)-(10*(max(mag)-min(mag))/18),'z-Sloan',fontproperties=font,fontsize = 12, color = 'brown')
text(0.3, max(mag)-(11*(max(mag)-min(mag))/18),'J-2MASS',fontproperties=font,fontsize = 12, color = '#6F4E37')
text(0.3, max(mag)-(12*(max(mag)-min(mag))/18),'H-2MASS',fontproperties=font,fontsize = 12, color = '#B87333')
text(0.3, max(mag)-(13*(max(mag)-min(mag))/18),'K-2MASS',fontproperties=font,fontsize = 12, color = '#827B60')
text(0.3, max(mag)-(14*(max(mag)-min(mag))/18),'uvw1-UVOT',fontproperties=font,fontsize = 12, color = '#7FFFD4')
text(0.3, max(mag)-(15*(max(mag)-min(mag))/18),'uvm2-UVOT',fontproperties=font,fontsize = 12, color = '#6960EC')
text(0.3, max(mag)-(16*(max(mag)-min(mag))/18),'uvw2-UVOT',fontproperties=font,fontsize = 12, color = '#7D0552')
text(0.3, max(mag)-(17*(max(mag)-min(mag))/18),'NUV-GALEX',fontproperties=font,fontsize = 12, color = '#95B9C7')
text(0.3, max(mag)-(18*(max(mag)-min(mag))/18),'FUV-GALEX',fontproperties=font,fontsize = 12, color = '#4B0082')
xlim(xl[0],xl[1])
ylim(yl[0],yl[1])

xlabel('Spectrum',size=18)
ylabel('Magnitude',size=18)
ax.minorticks_on()
show()
####################
##### writing legend on the file
####################
filekc.write("\n# ----------------------------------------------------------------------------------\n")
filekc.write("# Legend for SMS mode:\n")
filekc.write("# specTOspec        \t--> Magnitude computed with original spectrum\n")
filekc.write("# Hybrid     spec_BB\t--> Magnitude computed with a SN+Bbody hybrid\n")
filekc.write("# ----------------------------------------------------------------------------------\n")

sltime = _sleepc*len(snlist)
print '######################################################################'
print ''
print ' Evaluation done in %.0is, of which %.0is to take a nap to let rest your Random Access Memory (RAM) ' % (time,sltime)	
print ''
print ' \033[46mList of Magnitudes\033[0m ' , mag
print ' \033[44mVersion used\033[0m ' , method
print ''
print ' A text file has been created ===> Magnitudes_%.0i_%.0i.txt ' % (Tnowd,Tnow)		

