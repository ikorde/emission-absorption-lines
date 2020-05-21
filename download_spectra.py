import subprocess
import re

file = open(r'/Users/ikorde/Documents/caltech/deimos_redshift_linksIRSA.tbl.txt')

# find names of all spectra
addr_list = re.findall(r'<a href="/data/COSMOS/spectra/deimos/spec1d_fits/(.*?)">1D Fits spectrum</a>', file.read()) 

# download all spectra
for i in range(0,len(addr_list)):
	PATH = '/Users/ikorde/Documents/caltech/DEIMOS_spectra/' + addr_list[i]
	URL = 'https://irsa.ipac.caltech.edu/data/COSMOS/spectra/deimos/spec1d_fits/' + addr_list[i]
	cmd = 'wget -O %s %s ' % (PATH , URL)	
	print(cmd)
	subprocess.run(cmd, shell=True)
