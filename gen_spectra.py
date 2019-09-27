from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os,gc,pickle, re

genpath = '/Users/ikorde/Documents/caltech/'


""" Dataframe to store spectra from downloaded folder"""
filepathFits = genpath+'DEIMOS_spectra/spec1d.zd-1.092.hzoii-548248.fits'
ff = fits.open(filepathFits)
sd = ff[1].data
colnames = sd.columns.names
d = pd.DataFrame({'name':'spec1d.zd-1.092.hzoii-548248.fits', colnames[0]: [sd[0][0]], colnames[1]: [sd[0][1]], colnames[2]: [sd[0][2]] })

path = '/Users/ikorde/Documents/caltech/DEIMOS_spectra'
list_contents = os.listdir(path)

for i in range(1,len(list_contents)):
	sd = fits.open(path+'/'+list_contents[i])
	fd = sd[1].data
	colnames = fd.columns.names
	d = d.append({'name':list_contents[i], colnames[0]: [fd[0][0]] , colnames[1]: [fd[0][1]] , colnames[2]: [fd[0][2]] }, ignore_index=True)
	sd.close()

gc.collect()
d.to_pickle('deimos_spectra_new.pkl')


""" Dataframe to store IRSA table """
filepath = genpath+'deimos_redshift_linksIRSA.tbl.txt'
with open(filepath) as dradc_data: 
	lines = dradc_data.readlines()
lines = lines[4:]
ls = []
ls = [re.split('  +', x.strip()) for x in lines ]
df = pd.DataFrame(ls, columns = [ "id", "ra", "dec", "sel", "imag", "kmag", "zspec", "Qf", "Q", "remarks", "fits", "ascii", "jpeg", "2D", "temp"])
df.to_pickle('deimos_redshift_linksIRSA.pkl')


