import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from os.path import exists

from astropy import units as u
from astropy.io import fits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import medfilt

from specutils import Spectrum1D, SpectralRegion
from specutils.fitting import find_lines_derivative, find_lines_threshold, fit_lines, fit_generic_continuum
from specutils.manipulation import noise_region_uncertainty
from specutils.analysis import snr, equivalent_width


EAlines = {
	'LYA' 	: 1215.24,
	'HeII' 	: 1640.4,
	'OII' 	: 3727.092,	
	'HBeta' : 4862.68,
	'NII' 	: 6549.86,	
	'HAlpha': 6564.61,	
	'NII_2' : 6585.27, 
	'CIV'	: 1549,
	'OIII'	: 5008
}


def gauss(x,a,x0,sigma,c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

def error_estimate(flux,wvln,flux_sd,start,end,ycont): 
	fluxes_new = np.random.normal(loc=flux, scale=flux_sd, size=len(flux))
	return (np.trapz( fluxes_new[start:end] , x=wvln[start:end] ) / ycont) 


"""  Open DEIMOS Spectra """ 
deimos_spectra = pd.read_pickle('deimos_spectra.pkl')

""" Load DEIMOS RA/DC data from text to dataframe """ 
deimos_redshift_links = pd.read_pickle('deimos_redshift_linksIRSA.pkl')

""" Get specific item from deimos_redshift_linkIRSA """ 
sids = input("Enter list of spectra ID (ex. L665265, L539609) : ") #'L665265' #'L539609' #'L429443'
thres = input("Enter window for line: ")

breakpoint()
for sid in sids: 
	z = float(deimos_redshift_links.loc[deimos_redshift_links['id'] == sid, 'zspec'].item())
	n = deimos_redshift_links.loc[deimos_redshift_links['id'] == sid, 'fits'].item()
	RA = deimos_redshift_links.loc[deimos_redshift_links['id'] == sid, 'ra'].item()
	DEC = deimos_redshift_links.loc[deimos_redshift_links['id'] == sid, 'dec'].item()
	sname = n.replace('<a href="/data/COSMOS/spectra/deimos/spec1d_fits/', '')
	sname = sname.replace('">1D Fits spectrum</a>', '')

	""" Get flux and lambda from spectra file we chose to open """ 
	spec = sname #'spec1d.mn44.053.serendip.fits' #'spec1d.COS-23a.004.COSMOS_acal.fits' #'spec1d.COS-25.026.Rd-816509_acal.fits' #'spec1d.m25cc.044.civ_2862.fits' 
	flux = deimos_spectra.loc[ deimos_spectra['name'] == spec, 'fluxes'].tolist()[0][0]
	wvln = deimos_spectra.loc[ deimos_spectra['name'] == spec, 'lambdas'].tolist()[0][0]

	""" Convert data to Spectrum1D object for specutils methods """ 
	spectrum = Spectrum1D(flux*u.ct, spectral_axis=wvln*u.Angstrom)

	""" Lines using zero-crossings and first derivatives
		lines = find_lines_derivative(spectrum, threshold=5)

		Find lines by finding deviations larger than the noise_factor
		This method only works with continuum-subtracted spectra and the uncertainty must be defined on the spectrum
		noise_region_uncertainty does the above for us
	""" 
	noise_region = SpectralRegion(wvln[0]*u.Angstrom, wvln[1000]*u.Angstrom)
	spectrum = noise_region_uncertainty(spectrum, noise_region)

	""" Find lines using a noise threshold found above """ 
	lines = find_lines_threshold(spectrum, noise_factor = 1)
	lines = lines[0][:]

	""" Search for the E/A Lines we are looking for from the lines returned by specutils """ 
	found_lines = {}
	for line in lines: 
		for n,e in EAlines.items():
			if line.value > e*(1+z)-thres and line.value < e*(1+z)+thres:
				found_lines.update({ n: e*(1+z) })

	""" Continuum  """ 
	cont = fit_generic_continuum(spectrum)
	ycont = cont(wvln*u.Angstrom)
	f_interp_cont = interp1d(wvln, ycont)


	""" Create a mask for the lines that were found """ 
	f_interp_flux = interp1d(wvln, flux)

	mask_wvln = np.linspace(wvln[0], wvln[8191], 8192)
	mask_flux = ycont.copy() 	#[0]*8192


	""" Apply a median filter to the entire spectra """ 
	filtered_spec = medfilt(flux,3)*u.ct
	f_interp_filter = interp1d(wvln, filtered_spec)

	""" For each line that was found, label it on the plot, 
		Output the data,
		Find the index of the wvln where the line is centered,
		create a mask for +-30A next to that value from the filtered data, 
		Fit the Guassian to measure the lines
	""" 

	""" Plot the orginal spectra data """ 
	i=1
	flux_sd = np.std(flux[5700:7000])
	ILF = []

	if exists('FOUND_LINES.pkl'):
		df = pd.read_pickle('FOUND_LINES.pkl')
	else :
		df = pd.DataFrame(columns=EAlines.keys(), index=[sid])


	df_curr = pd.DataFrame(index=["z-spec", "RA", "DEC", "Flux Lower", "Flux Mean", "Flux Upper", "EW Lower", "EW", "EW Upper"],columns=EAlines.keys())
	df_curr.at['z-spec'] = z

	for n,e in found_lines.items():
		""" search for where the wavelength=center of line (e)"""
		ind = np.where((wvln <= int(e)+1) & (wvln >= int(e)))[0][0]

		"""window +-40A for dist of center of line"""
		start = ind-40 
		end	  = ind+40
		inds = np.linspace(e-40, e+40, 1000)

		mask_flux[start:end] =filtered_spec[start:end] #f_interp_filter(inds)

		""" chopped data for fiding fit """
		chop_flux = mask_flux[start:end]
		chop_wvln = mask_wvln[start:end] 
		
		"""convert data to dimensionless np array for fitting guassian"""
		chop_flux = np.array(chop_flux*u.dimensionless_unscaled)
		chop_wvln = np.array(chop_wvln*u.dimensionless_unscaled)

		"""fit guassian using scipy curve fit"""
		LINECENTER = mask_wvln[ind]
		LINEFWHM = 5*2.35
		popt , pcov = curve_fit(gauss,xdata=chop_wvln,ydata=chop_flux,
			p0=[np.nanmax(chop_flux) , LINECENTER , LINEFWHM/2.35 , 1],
			bounds=(np.array([-np.inf, LINECENTER-2,0,-2]) , 
			np.array([np.inf, LINECENTER+2, np.inf, 1]) ), )

		mu_fit = popt[1]
		sd_fit = popt[2]

		Y_fit = gauss(chop_wvln , *popt)
		f_interp_yfit = interp1d(chop_wvln,Y_fit)

		"""EW"""
		xlower = mu_fit - 0.41*sd_fit
		xupper = mu_fit + 0.41*sd_fit

		EW_lower = np.trapz( Y_fit, x=chop_wvln )/f_interp_cont(xlower)
		EW_upper = np.trapz( Y_fit, x=chop_wvln )/f_interp_cont(xupper)
		EW_center = np.trapz( Y_fit, x=chop_wvln )/f_interp_cont(e)

		"""" Estimate Error 
			 Generate 1000 random fluxes based on SD of original data 
			 Find the ILF and EW of each of these trials for each line in the spectrum 
		"""
		ILF.append(sorted(list ([ error_estimate(flux,wvln,flux_sd,start,end, np.array(f_interp_cont(e)*u.dimensionless_unscaled)) for j in range(0,1000)] ) ))

		""" Update the dataframes """
		df.at[sid, n] = 'X'
		df_curr.at['Flux Upper', n] = xupper
		df_curr.at['Flux Mean', n] = mu_fit
		df_curr.at['Flux Lower', n] = xlower
		df_curr.at['EW', n] = EW_center
		df_curr.at['EW Lower', n] = EW_lower
		df_curr.at['EW Upper', n] = EW_upper
		df_curr.at['RA', n] = RA
		df_curr.at['DEC', n] = DEC

		i+=1

	df_curr.to_excel(r'/Users/ikorde/Documents/caltech/data/%s_data.xlsx'%(sid))


df.to_excel(r'/Users/ikorde/Documents/caltech/data/all_spectralines.xlsx')

# df.to_pickle('FOUND_LINES.pkl')


