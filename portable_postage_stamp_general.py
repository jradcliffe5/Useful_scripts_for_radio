from astropy.table import Table
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import pickle
import pandas as pd
import os, sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
from astropy.coordinates import SkyCoord
from matplotlib import rcParams
import traceback
import matplotlib.gridspec as gridspec


def cutout_sources_general(fitsfile,ra_col,dec_col,col_unit,catalogue,subimsize, logthresh,log_contour_scale,nlevs,flux_scaler,subplot_x,subplot_y,plot_prefix,z_col):
	matplotlib.rcParams.update({'font.size': 10})
	logthresh = -1*logthresh
	nlevs = nlevs
	do_rms_contours = 'False'
	subimsize = int(subimsize / 2.)
	hdu = fits.open(fitsfile)
	image_data = hdu['PRIMARY'].data
	try:
		bmaj = hdu[0].header['BMAJ'] / hdu[0].header['CDELT2']
		bmin = hdu[0].header['BMIN'] / hdu[0].header['CDELT2']
		bpa = hdu[0].header['BPA']
		beam_info = True
	except KeyError:
		beam_info = False
	ndims = image_data.ndim
	if ndims > 2:
		image_data = image_data[0, 0, :, :]
	c = SkyCoord(catalogue[ra_col],catalogue[dec_col],unit=col_unit)
	print('Making postage stamps for the %d sources in the catalogue' % len(catalogue[ra_col]))
	page_counter = 0
	running_counter = 0
	for i in range(len(catalogue[ra_col])):
		#print(running_counter % (subplot_x*subplot_y))
		if running_counter % (subplot_x*subplot_y) == 0:
			plt.clf()
			fig = plt.figure(figsize=(4*subplot_x, 4*subplot_y))
			gs = gridspec.GridSpec(subplot_x, subplot_y)
			gs.update(wspace=0.11, hspace=0.15)
			x_counter = 0
			y_counter = 0
			page_counter +=1
		try:
			wcs = WCS(hdu['PRIMARY'].header)
			if ndims > 2:
				central_pix_coord = wcs.wcs_world2pix(c.ra[i],c.dec[i],1,1,1)
			else:
				central_pix_coord = wcs.wcs_world2pix(c.ra[i],c.dec[i],1)
			'''
			if units == 'uJy':
				flux_scaler = 1E6
			elif units == 'mJy':
				flux_scaler = 1E3
			elif units == 'Jy':
				flux_scaler = 1
			else:
				flux_scaler = flux_scaler
			'''
			flux_scaler = flux_scaler
			RA_min2 = int(central_pix_coord[0]) - subimsize
			RA_max2 = int(central_pix_coord[0]) + subimsize
			Dec_min2 = int(central_pix_coord[1]) - subimsize
			Dec_max2 = int(central_pix_coord[1]) + subimsize
			### Cut out image
			image_data2 = image_data[Dec_min2:Dec_max2,RA_min2:RA_max2] * flux_scaler
			RA_pix = (RA_min2 + RA_max2) / 2.
			Dec_pix = (Dec_min2 + Dec_max2) / 2.
			## Adjust wcs
			wcs2 = wcs
			#print(x_counter,y_counter)
			print(x_counter,y_counter)
			ax = fig.add_subplot(gs[x_counter,y_counter],projection=wcs2)
			if ndims>2:
				RA_w, Dec_w = wcs.wcs_pix2world(RA_pix, Dec_pix, 1,1,1)[0:2]
				wcs2.wcs.crval = [RA_w, Dec_w,wcs.wcs.crval[2],wcs.wcs.crval[3]]
				wcs2.wcs.crpix = [subimsize, subimsize,wcs.wcs.crpix[2],wcs.wcs.crpix[3]]
			else:
				RA_w, Dec_w = wcs.wcs_pix2world(RA_pix, Dec_pix, 1)
				wcs2.wcs.crval = [RA_w, Dec_w]
				wcs2.wcs.crpix = [subimsize, subimsize]
			### Make name
			coord = SkyCoord(
				catalogue[ra_col].iloc[i], catalogue[dec_col].iloc[i], unit=col_unit)
			if coord.dec.dms.d < 0:
				neg = '-'
			else:
				neg = '+'
			if len(str(coord.ra.hms.s).split('.')[0]) == 1:
				ra_s = '0%2.2f' % coord.ra.hms.s
			else:
				ra_s = '%2.2f' % coord.ra.hms.s
			if len(str(coord.dec.dms.s).split('.')[0]) == 1:
				dec_s = '0%2.2f' % coord.dec.dms.s
			else:
				dec_s = '%2.2f' % coord.dec.dms.s
			name = 'J%s%s%s%s%s%s%s' % ('%02d' % int(coord.ra.hms.h),
										'%02d' % int(coord.ra.hms.m),
										'%s' % ra_s,
										neg, '%02d' % int(coord.dec.dms.d),
										'%02d' % int(coord.dec.dms.m),
										'%s' % dec_s)
			print('%d) Plotting %s' % (i+1,name))
			#ax = plt.subplot(projection=wcs2)
			### Set coordinate formats
			lon = ax.coords['ra']
			lat = ax.coords['dec']
			lon.set_major_formatter('d.ddd')
			lat.set_major_formatter('d.ddd')
			'''
			if x_counter==0:
				lon.set_major_formatter('d.ddd')
			else:
				lon.set_ticklabel_visible(False)
				lon.set_axislabel('')
			### Some VLBI specific shit
			if (y_counter + 1) == subplot_y:
				lat.set_major_formatter('d.ddd')
			else:
				lat.set_ticklabel_visible(False)
				lat.set_axislabel('')
			'''
			### Set labels for axes
			#lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
			#lat.set_axislabel('Declination (J2000)', minpad=1)
			lon.set_ticks(number=3)
			#ax.set_xlim(int(RA_pix[i]) - subimsize,int(RA_pix[i]) + subimsize)
			#ax.set_ylim(int(Dec_pix[i]) - subimsize,int(Dec_pix[i]) + subimsize)
			### Makes colorbar on top
			divider = make_axes_locatable(ax)
			cax = divider.append_axes(
				"top", size="5%", pad=0.00, axes_class=matplotlib.axes.Axes)
			if np.max(image_data2) > 60:
				im = ax.imshow(
					image_data2,
					origin='lower',
					cmap="magma",
					interpolation="bicubic",
					norm=matplotlib.colors.SymLogNorm(10**-logthresh))
				tick = [np.min(image_data2)]
				if log_contour_scale == 'e':
					tick = np.append(tick,
									 np.around(
										 np.logspace(
											 np.log(np.std(image_data2)),
											 np.log(np.max(image_data2)),
											 nlevs,
											 base=np.e,
											 endpoint=True),
										 decimals=0)[1:])
				elif log_contour_scale == '2':
					tick = np.append(tick,
									 np.around(
										 np.logspace(
											 np.log2(np.std(image_data2)),
											 np.log2(np.max(image_data2)),
											 nlevs,
											 base=2,
											 endpoint=True),
										 decimals=0)[1:])
				elif log_contour_scale == '10':
					tick = np.append(tick,
									 np.around(
										 np.logspace(
											 np.log10(np.std(image_data2)),
											 np.log10(np.max(image_data2)),
											 nlevs,
											 base=10,
											 endpoint=True),
										 decimals=0)[1:])
				else:
					print('log_contour_scale can only be \'e\', 2 or 10')
					exit()
				tick = tick.astype(int)
				cb = plt.colorbar(
					orientation="horizontal", mappable=im, cax=cax, ticks=tick)
			else:
				im = ax.imshow(
					image_data2,
					origin='lower',
					cmap="magma",
					interpolation="bicubic")
				tick = [np.min(image_data2)]
				tick = np.append(tick,
								 np.around(
									 np.linspace(
										 np.std(image_data2), np.max(image_data2),
										 nlevs),
									 decimals=1)[1:].astype(int))
				tick = np.append(tick, np.max(image_data2))
				tick = tick.astype(int)
				cb = plt.colorbar(
					orientation="horizontal",
					mappable=im,
					cax=cax,
					ticks=tick,
					format='{:.0f}')
			### Set tick position
			cb.ax.xaxis.set_ticks_position('top')
			'''
			if np.max(image_data2) > 60:
				if do_rms_contours == 'False':
					levs = []
				else:
					levs = [-1 * np.std(image_data2), np.std(image_data2)]
				if log_contour_scale == 'e':
					levs = np.append(levs,
									 np.around(
										 np.logspace(
											 np.log(np.std(image_data2)),
											 np.log(np.max(image_data2)),
											 nlevs,
											 base=np.e,
											 endpoint=True),
										 decimals=0)[1:-1])
				elif log_contour_scale == '2':
					levs = np.append(levs,
									 np.around(
										 np.logspace(
											 np.log2(np.std(image_data2)),
											 np.log2(np.max(image_data2)),
											 nlevs,
											 base=2,
											 endpoint=True),
										 decimals=0)[1:-1])
				elif log_contour_scale == '10':
					levs = np.append(levs,
									 np.around(
										 np.logspace(
											 np.log10(np.std(image_data2)),
											 np.log10(np.max(image_data2)),
											 nlevs,
											 base=10,
											 endpoint=True),
										 decimals=0)[1:-1])
				else:
					print('log_contour_scale can only be \'e\', 2 or 10')
					exit()
			else:
				levs = [-1 * np.std(image_data2), np.std(image_data2)]
				levs = np.append(levs,
								 np.around(
									 np.linspace(
										 np.std(image_data2), np.max(image_data2),
										 nlevs),
									 decimals=1)[1:-1])
			'''
			#cont = ax.contour(image_data2, levels=levs, cmap='gray_r', alpha=0.5)
			ax.scatter(coord.ra.degree, coord.dec.degree, transform=ax.get_transform('icrs'), s=60, edgecolor='k', facecolor='k', marker='+', alpha=1)
			#cb.add_lines(cont)
			cb.ax.set_xticklabels(tick)
			cb.ax.minorticks_off()
			if beam_info == True:
				circ = Ellipse(
					(10, 10),
					width=bmin,
					height=bmaj,
					angle=bpa,
					fill=False,
					color='w',
					hatch='xxxxxx')
				ax.add_patch(circ)
			ax.set_zorder(20)
			ax.text(
				0.1,
				0.9,
				r'{\bf %s}' % name,
				verticalalignment='center',
				horizontalalignment='left',
				transform=ax.transAxes,color='w')
			if z_col !=False:
				ax.text(
					0.1,
					0.85,
					r'{\bf $z=%.2f$}' % catalogue[z_col].iloc[i],
					verticalalignment='center',
					horizontalalignment='left',
					transform=ax.transAxes,color='w')
			if ((subplot_x==(x_counter+1)) and (subplot_y==(y_counter+1))) or (i==(len(catalogue[ra_col])-1)):
				plt.savefig('%s_%d.pdf' % (plot_prefix,page_counter),dpi=100, bbox_inches='tight', clobber=True,orientation='landscape')
				plt.clf()
			if (int(x_counter+1) % subplot_x) == 0:
				y_counter += 1
				x_counter = 0
			else:
				x_counter += 1
			del wcs2
			running_counter += 1
		except ValueError:
			traceback.print_exc()
	hdu.close()
	del image_data
	os.system('tar cf %s.tar.gz %s_*.pdf' % (plot_prefix,plot_prefix))
	#os.system('rm postage_stamps/%s.pdf' % fitsfile.split('.fits')[0])
	#os.system('ulimit -S -n 16192')
	#os.system('pdfjoin --no-landscape --rotateoversize False -o %s.pdf *pdf' % fitsfile.split('.fits')[0])
	#os.system('gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r50 -sOutputFile=%s_small.pdf %s.pdf' % ( fitsfile.split('.fits')[0],fitsfile.split('.fits')[0]))
	#os.system('mv *.pdf postage_stamps/')
