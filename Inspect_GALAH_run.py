from Inspect_GALAH_class import *
from glob import glob

sobjects = list([])
spectra_dir = 'SPECTRA/dr5.2/'
for dir in glob(spectra_dir+'*'):
	date = (dir.split('/')[-1])
	for fit in glob(spectra_dir+date+'/standard/com/*.fits'):
		sobjects.append(fit.split('/')[-1][:-6])

# Get unique s_object ids
sobjects = np.unique(sobjects)

for s_id in sobjects:
	# Now let's create the class FITS for a given sobject_id and use the provided functions on it!
	fits = fits_class(sobject_id=s_id)

	# ## Assuming you want to look at specific lines (and check RV), then use define the 'help_lines' variable before plotting, otherwise define 'help_line=False'
	help_lines=False
	help_lines=[
	    (r'$\mathrm{H_\alpha}$',6562.7970),
	    (r'$\mathrm{H_\beta}$' ,4861.3230),
	    (r'$\mathrm{Li}$'      ,6707.7635)
	    ]


	# ## Plot normalised spectrum on 4 axes
	fits.plot_norm_spectrum_on4axes(help_lines=help_lines)

	# ## Plot raw spectrum on 4 axes
	fits.plot_raw_spectrum_on4axes(help_lines=help_lines)

	# ## Plot Hbeta / Halpha / Li windows
	fits.plot_HbHaLi()

	# # Routine to get SME data in class 'sme'
	sme = sme_class(sobject_id=s_id)

	for each_mode in sme.sme['mode']:
	    sme.plot_mode(each_mode)
