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


def gauss(x,a,x0,sigma,c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

def gauss_fit(wvln, flux, linecenter, linestd):
	LINECENTER = linecenter 	
	LINEFWHM = linestd #5*2.35
	popt , pcov = curve_fit(gauss,xdata=chop_wvln,ydata=chop_flux,
		p0=[np.nanmax(chop_flux) , LINECENTER , LINEFWHM/2.35 , 1],
		bounds=(np.array([-np.inf, LINECENTER-10,0,-2]) , 
		np.array([np.inf, LINECENTER+10, np.inf, 1]) ), )

	mu_fit = popt[1]
	sd_fit = popt[2]

	Y_fit = gauss(chop_wvln , *popt)
	return(Y_fit, mu_fit, sd_fit)

def error_estimate(flux,wvln,flux_sd_e,center,start,end,ycont): 
	fluxes_new = (np.random.normal(loc=flux, scale=np.abs(flux_sd_e), size=len(flux)))
	f_interp_fluxes_e = interp1d(wvln, fluxes_new)
	Y_fit_e, mu_fit_e, sd_fit_e = gauss_fit(chop_wvln, fluxes_new[start:end], center, np.std(chop_wvln))
	Y_fit_e = np.where(Y_fit_e < ycont(mu_fit_e), ycont(mu_fit_e), Y_fit_e)
	if mu_fit_e > max(wvln): return 
	flux_center_e = float(f_interp_fluxes_e(mu_fit_e))
	EW_center_e = np.trapz( Y_fit_e, x=chop_wvln )/ycont(mu_fit_e)
	return [Y_fit_e,mu_fit_e, flux_center_e, EW_center_e]



""" Initializations """
all_spectra = False
show_graphs = True
num_error_trials = 500
cropped_window = 20
line_window = 30
EAlines = {
	'LYA' 	: 1215.24,
	'HeII' 	: 1640.4,
	'OII' 	: 3727.092,	
	'HBeta' : 4862.68,
	'NII' 	: 6549.86,	
	'HAlpha': 6564.61,	
	'CIV'	: 1549,
	'OIII'	: 5008
}

"""  Open DEIMOS Spectra """ 
deimos_spectra = pd.read_pickle('deimos_spectra.pkl')

""" Load DEIMOS RA/DC data from text to dataframe """ 
deimos_redshift_links = pd.read_pickle('deimos_redshift_linksIRSA.pkl')

""" Get specific item from deimos_redshift_linkIRSA """ 
if all_spectra: 
	sids = deimos_redshift_links.loc[:, 'ID'].tolist()

else: 
	sids = input("Enter list of spectra ID (ex. L665265, L539609) : ") #'L665265' #'L539609' #'L429443'
	sids = sids.split(', ')
	print('\n')

if exists('df_all.pkl'):
	df_all = pd.read_pickle('df_all.pkl')
else :
	lines = np.repeat([*EAlines.keys()],10)
	names = ['SNR', 'linelow', 'linecenter', 'lineupper', 'fluxlow', 'fluxcenter', 'fluxupper', 'EWlow', 'EWcenter', 'EWupper']*len(EAlines.keys())
	index = pd.MultiIndex.from_arrays([lines,names])
	df_all = pd.DataFrame(columns=index)	#pd.DataFrame(columns=['Z', 'RA', 'DEC', *EAlines.keys()])
	df_all.insert(0,'Z',True)
	df_all.insert(1,'RA',True)
	df_all.insert(2,'DEC',True)
	# df_all.set_index('Name', inplace=True)

for sid in sids:
	# try:
	if pd.isna(sid): 
		sid = 'N/A'

	try: 
		z = float(deimos_redshift_links.loc[deimos_redshift_links['ID'] == sid, 'zspec'].item())
	except ValueError: 
		print('Invalid redshift for galaxy %s, continuing... \n'%sid)
		continue
	n = deimos_redshift_links.loc[deimos_redshift_links['ID'] == sid, 'fits1d'].item()
	if pd.isna(n): 
		print('No spectra for galaxy %s \n'%sid)
		continue

	RA = deimos_redshift_links.loc[deimos_redshift_links['ID'] == sid, 'ra'].item()
	DEC = deimos_redshift_links.loc[deimos_redshift_links['ID'] == sid, 'dec'].item()
	sname = n.replace('<a href="/data/COSMOS/spectra/deimos/spec1d_fits/', '')
	sname = sname.replace('">1D Fits spectrum</a>', '')
	
	""" Get flux and lambda from spectra file we chose to open """ 
	spec = sname #'spec1d.mn44.053.serendip.fits' #'spec1d.COS-23a.004.COSMOS_acal.fits' #'spec1d.COS-25.026.Rd-816509_acal.fits' #'spec1d.m25cc.044.civ_2862.fits' 
	flux = deimos_spectra.loc[ deimos_spectra['name'] == spec, 'FLUX'].tolist()[0][0]
	wvln = deimos_spectra.loc[ deimos_spectra['name'] == spec, 'LAMBDA'].tolist()[0][0]

	""" Convert data to Spectrum1D object for specutils methods 
		First apply median filter to flux
	""" 
	filtered_flux = medfilt(flux,3)
	spectrum = Spectrum1D(filtered_flux*u.ct, spectral_axis=wvln*u.Angstrom)

	""" Lines using zero-crossings and first derivatives
		lines = find_lines_derivative(spectrum, threshold=5)

		Find lines by finding deviations larger than the noise_factor
		This method only works with continuum-subtracted spectra and the uncertainty must be defined on the spectrum
		noise_region_uncertainty does the above for us
	""" 
	noise_region = SpectralRegion(wvln[1000]*u.Angstrom, wvln[-1000]*u.Angstrom)
	try: 
		spectrum = noise_region_uncertainty(spectrum, noise_region)
	except ValueError: 
		continue
	lines = find_lines_threshold(spectrum, noise_factor = 1)[0][:]

	""" Search for the E/A Lines we are looking for from the lines returned by specutils """ 
	found_lines = {}
	for line in lines: 
		for n,e in EAlines.items():
			if line.value > e*(1+z)-line_window and line.value < e*(1+z)+line_window :
				found_lines.update( {n: e*(1+z)} )

	df_all.loc[sid, 'Z'] = z
	df_all.loc[sid, 'RA'] = RA
	df_all.loc[sid, 'DEC'] = DEC

	""" If no lines found, continue """
	if not found_lines: 
		print('No lines in galaxy %s \n'%sid)
		continue

	""" Continuum  """ 
	cont = fit_generic_continuum(spectrum)
	ycont = cont(wvln*u.Angstrom)
	f_interp_cont = interp1d(wvln, ycont)

	""" Create a mask for the lines that were found """ 
	f_interp_flux = interp1d(wvln, flux)
	f_interp_filter = interp1d(wvln, filtered_flux)

	mask_wvln = wvln.copy() 	
	mask_flux = ycont.copy() 	

	""" Plot the orginal spectra data """ 
	if(show_graphs):
		plt.figure(0)
		plt.plot(wvln, flux, 'lightgray')
		plt.plot(wvln, filtered_flux)
	i=1
	flux_sd = np.std(filtered_flux)
	all_data = []

	""" For each line that was found, create a mask for +-windowA next to that value from the filtered data, 
		fit the Guassian to measure the lines, calculate error on fit
	""" 
	# print('E/A lines in galaxy %s:'%sid)
	for n,el in found_lines.items():
		# print(n,'at wavelength', el)

		""" search for where the wavelength=center of line (e)"""
		ind = np.where((wvln <= int(el)+1) & (wvln >= int(el)))[0][0]

		"""window for dist of center of line"""
		start = ind - cropped_window 
		end	  = ind + cropped_window

		mask_flux[start:end] =filtered_flux[start:end]*u.ct #f_interp_filter(inds)

		""" chopped data for fiding fit """
		chop_flux = mask_flux[start:end]
		chop_wvln = mask_wvln[start:end] 
		
		"""convert data to dimensionless np array for fitting guassian"""
		chop_flux = np.array(chop_flux*u.dimensionless_unscaled)
		chop_wvln = np.array(chop_wvln*u.dimensionless_unscaled)

		"""Dictionary to store analysis values and load to datatable"""
		# keys = ['Flux Lower', 'Flux', 'Flux Upper','Line Lower', 'Line', 'Line Upper', 'EW Lower', 'EW', 'EW Upper']
		# line_ew_flux = {}
		# line_ew_flux.fromkeys(keys,None)

		"""fit guassian using scipy curve fit"""
		try: 
			Y_fit, mu_fit, sd_fit = gauss_fit(chop_wvln, chop_flux, el, np.std(chop_wvln))
		except RuntimeError: 
			continue
		f_interp_yfit = interp1d(chop_wvln,Y_fit)

		"""EW"""
		# xlower = mu_fit - 0.41*sd_fit
		# xupper = mu_fit + 0.41*sd_fit

		# breakpoint()
		# EW_lower = np.trapz( Y_fit, x=chop_wvln )/f_interp_cont(xlower)
		# EW_upper = np.trapz( Y_fit, x=chop_wvln )/f_interp_cont(xupper)
		EW_center = np.trapz( Y_fit, x=chop_wvln )/f_interp_cont(el)
		"""" Estimate Error 
			 Generate 1000 random fluxes based on SD of original data 
			 Find the ILF and EW of each of these trials for each line in the spectrum 
		"""
		
		all_data.append(list(error_estimate(filtered_flux, wvln, flux_sd, el, start, end, f_interp_cont) for j in range(0,num_error_trials)))
		all_data = all_data[0]
		if not any(all_data): continue 

		tmp = np.percentile(list(zip(*all_data))[1],(50,50-68/2,50+68/2))
		linelow = tmp[0]
		linecenter = tmp[1]
		lineupper = tmp[2]

		tmp1 = np.percentile(list(zip(*all_data))[2],(50,50-68/2,50+68/2))
		fluxlow = tmp1[0]
		fluxcenter =  tmp1[1]
		fluxupper = tmp1[2]

		tmp2 = np.percentile(list(zip(*all_data))[3],(50,50-68/2,50+68/2))
		EWlow = tmp2[0]
		EWcenter = tmp2[1]
		EWupper = tmp2[2]

		SNR = fluxcenter/flux_sd

		l = [SNR, linelow, linecenter, lineupper, fluxlow,fluxcenter,fluxupper,EWlow,EWcenter,EWupper]
		df_all.at[sid, n] = l

		y_fit_errorfix = list(zip(*all_data))[0]

		""" Plot data """
		if(show_graphs):
			plt.figure(0)
			plt.annotate(n, (el, 10), ha='left')
			fig, (fit,error) = plt.subplots(2)
			plt.figure(i)
			fit.plot(chop_wvln, flux[start:end], 'lightgray', label='original')
			fit.plot(chop_wvln, chop_flux, label='error fixed')
			fit.plot(chop_wvln, y_fit_errorfix[0], 'orange', label='fit') #gauss(X_fit, *popt))
			# fit.plot(linelow, f_interp_yfit(linelow), 'r*')
			# fit.plot(lineupper, f_interp_yfit(lineupper), 'r*')
			# fit.axvspan(linecenter-(EWcenter-(EWcenter/2)), linecenter+(EWcenter+(EWcenter/2)), alpha=0.5)
			fit.legend()
			fit.title.set_text('%s: %s Fit'%(sid,n))
			# error.hist(list(zip(*all_data))[2], range=[EWlow,EWupper], bins=10)
			# error.title.set_text('Distribution of equivalent widths over %i random trials'%num_error_trials)
			# plt.tight_layout()
		all_data = []
		i+=1
	# except Exception as e:
	# 	print(e)
	# 	breakpoint()

	print('')
	""" Plot the remaining data """ 
	if(show_graphs):
		plt.figure(0)
		plt.plot(mask_wvln,	mask_flux, color="orange")
		plt.xlabel( 'Wavelength (Angstrom)' )
		plt.ylabel( 'Flux (Counts)' )
		plt.title(sid)
		# plt.plot(wvln, ycont, color="green")
		plt.show()

df_all.to_pickle('df_all.pkl')
df_all.to_excel('df_all.xlsx')

print('Data saved.')

breakpoint()




