###########################################################################

							  S3 package
							 ------------
							  S.N.A.K.E.
							   S.N.A.P.
							    S.M.S       

					C. Inserra  v1.1.0 29/10/2015

###########################################################################
The package contains:
- STHREE (general descriptor)
- SNAKE
- SNAKELOOP
- SNAP
- SMS

Usage:
> SNAKE (in terminal wherever you want)

Help:
SNAKE -h
SNAKELOOP -h
SNAP -h
SMS -h

The help function shows the correct usage in the terminal and the other
options available such as change the error on the redshift (-r) or the
sleeping time between loops (-s) for those codes working with loops.
###########################################################################

----------------------------- S.N.A.K.E. ----------------------------------

The SuperNova Algorithm for K-correction Evaluation (S.N.A.K.E.) 
is based on the formula m(x) = M(y) + DM + K(y,x),
then if you want to add the correction to the observed magnitude you 
have to flip the sign                         

The user can introduce the values or hit enter. In the last case the
values between brackets will be taken as default value.

The SNAKE_loop script works in the same way of the SNAKE one, but it will 
plot the values and create a text file with all the relevant information.

Errors are now introduced and are treated as a r.m.s of four different
errors considered individually as a direct measurements errors, hence no
monte-carlo simulation to create a statistical error has been executed.
The error is due to five terms (see Inserra & Smartt 2014 or Inserra et 
al. 2015). Here the errors regarding a template spectrum are not included,
while the errors on the redshift, the fluxes zero points for the integration
of the filter over the spectrum range are accounted for. This together 
with an additional error due to the assumption that your spectrum is similar 
to a black-body function in the case part of your filter does not cover the
real spectrum.  

BE SURE that the spectrum is at least calibrated in flux relatively (e.g.
the synthetic colours match those from photometry).

If the spectrum is also calibrated in flux absolutely, the magnitudes 
values retrieved by the script should be correct and in agreement with 
those from photometry.

------------------------------ S.N.A.P. -----------------------------------

The SuperNova Algorithm for P-correction (S.N.A.P) apply a passband
correction that is the "young sibling" of the more famous S-correction
(Stritzinger et. al 2002, Pignata et al. 2003). It follows the simple
mathematical formula:

P(lambda) = F(lambda) * QE(lambda)

Where F(lambda) is the filter function and QE(lambda) is the quantum 
efficiency of the telescope used. Note that the original S-correction
had three terms. The first one is related to the continuum atmospheric
transmission profile of the site that it is not included in the 
P-correction since SN magnitudes are evaluated through sequence stars 
calibrated with Sloan stars (and hence such correction is already taken 
in account at this level) or with Landolt stars and subsequent use of 
programmes that apply such correction. If none of these two method is used 
you can find on-line such extinction curves for several telescopes. 
The second term is related to the mirror reflectivity function that is 
usually around 90% with a small depression around 8000 Angs. Such function
can change from telescope to telescope but these changes are less than 3%
and they affects the final measurements of the order of the photometric
errors, thus - given the difficulty to retrieve it for all the telescope -
the P-correction does not take it in account. The last term is the lens 
throughput that is almost constant through the optical regime for each 
telescope and hence it is not taken in account. This might cause some 
problems at wavelengths bluer than 3400 Angs, i.e. affecting the U passband.

The usage is based on the formula m(P) = m(F) - P(lambda), where m(P)
is the magnitude of the instrumental passband and m(F) is the 
magnitude you already evaluated from your photometry at Telescope XXX 
(XXX = NTT or NOT or TNG or PS1 or LT or LSQ or ASIAGO or SKYMAPPER 
or LCOGT or OGLE)

The script will plot the values and create a text file with all the 
relevant information.

Errors are now included and are evaluated like in SNAKE. Here we do not
have an error on the redshift.

BE SURE that the spectrum is at least calibrated in flux relatively (e.g.
the synthetic colours match those from photometry).

If the spectrum is also calibrated in flux absolutely, the magnitudes 
values retrieved by the script should be correct and in agreement with 
those from photometry.

------------------------------ S.M.S. -------------------------------------

The Synthetic Magnitudes from Spectra (S.M.S.) code is a python version,
IRAF-free, of the STSDAS/HST tool calcphot. It works with .fits and it
will plot the values and create a text file with all the relevant 
information.

###########################################################################

					ACKNOWLEDGMENTS AND FAQ

Many thanks to Elisabeth E. E. Gall for testing and debugging.  

I am also indebted to David R. Young and Ken W. Smith for their assistance 
with Python.

Please report any problems to c.inserra@qub.ac.uk  
or open a ticket on https://github.com/cinserra (best way)

###########################################################################

						RELEASE NOTES

Please, read the ReleaseNotes.txt file if you want to know what has been
changed/updated/removed/fixed from the previous version(s)

###########################################################################

						INSTALLATION

S3 is written in python and requires the following package:

- Python 2.5 or Python 2.6 or Python 2.7
   these modules have to be installed:
        - numpy
        - scipy
	    - pyraf
	    - matplotlib
	    - pyfits
- Iraf 

If you have UREKA installed is even better. Otherwise you can
retrieve it here: http://ssb.stsci.edu/ureka/

WARNING:
- Occasionally could happen that the computer is not able to properly
call the programmes from your bin. Check where they are (which SNAKE)
and be sure that the first line is #!/usr/bin/env python	

###########################################################################
1) extract the files from the tarball
> tar -xvf S3.tar
> cd S3

2) install the programs as user

YOU DO NOT NEED TO BE ROOT!!!!
Between parentheses there are optionals commands

> python setup.py install  (--prefix=/home/user/s3) (--record file.txt)

It is preferable to install it under your user directory (in a OSX
system it will be /Users/nameoftheuser/s3)

You can also install the scripts in the same directory where you downloaded
the file avoiding the --prefix=/home/user/xxx 

file.txt will be written, creating a log of the installation

###########################################################################

To remove a previous version 

3) As a first attempt try to install it again, it should work,
if the version doesn't match proceed to next step

4)
- delete the S3 directory in your site-package path
- delete the s3****.egg-info from the same directory
- delete all the executable: 
STHREE, SNAKE, SNAKELOOP, SNAP, SMS

5) to do so type the following commands:

> which SNAKE
(e.g. /Users/Cosmo/Ureka/variants/common/bin/SNAKE )
> rm -rf /Users/Cosmo/Ureka/variants/common/bin/SN* 
(check if you do not have other bin with SN*)
> rm -rf /Users/Cosmo/Ureka/variants/common/bin/SMS
> rm -rf /Users/Cosmo/Ureka/variants/common/bin/STHREE

6) open python:

> python

>>> import s3
>>> s3.__path__
(e.g. ['/Users/Cosmo/Ureka/variants/common/lib/python2.7/site-packages/s3'])

7) then on the terminal (or X11):

> rm -rf /Users/Cosmo/Ureka/variants/common/lib/python2.7/site-packages/s3
> rm -rf /Users/Cosmo/Ureka/variants/common/lib/python2.7/site-packages/s3*.egg-info

ALTERNATIVE:
If you used the option " --record files.txt " during the installation
you can run the following command in the terminal:
> cat files.txt | xargs sudo rm -rf


8) Now you can repeat points 1 to 2

###########################################################################
