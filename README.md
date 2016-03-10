# S3
Python suite containing four algorithms related to broadband filter observations of distant objects and/or observations with different set of filters

# SNAKE
he SuperNova Algorithm for K-correction Evaluation (S.N.A.K.E.) 
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

# SNAP
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

# SMS
The Synthetic Magnitudes from Spectra (S.M.S.) code is a python version,
IRAF-free, of the STSDAS/HST tool calcphot. It works with .fits and it
will plot the values and create a text file with all the relevant 
information.
