
#!/usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import datetime
import time
from sys import exit
from matplotlib import colors, pyplot as plt
from functools import reduce
import matplotlib.cm as cm
import seaborn as sns
from astropy.io import ascii, fits
from astropy.wcs import wcs
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.interpolate import interp2d
import matplotlib.mlab as mlab
import scipy, pylab
import rpy2
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import math
from pysextractor import SExtractor

__author__ = 'pnovais'

ini=time.time()

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#definindo a classe que ira ler as imagens fits
def get_image(f_sdss):
    img = f_sdss[0].data
#    sky = f_sdss[2].data
    return img

#abertura do arquivo com o nome das imagens, nas n bandas
df_fit = pd.read_csv('data/arquivo_fits.csv')

'''
================================================================================
Rodando o SExtractor na imagem na banda r, criando uma segmentacao e um catalogo
com os objetos obtidos
ATUALIZAR NOME DA BANDA DE SEGMENTACAO
================================================================================
'''
fname = 'data/frame-r-002507-4-0226.fits'
sex = SExtractor()
sex.config['PARAMETERS_LIST'].append('FLUX_ISO')
sex.config['PARAMETERS_LIST'].append('MAG_ISOCOR')
sex.config['PARAMETERS_LIST'].append('MAG_AUTO')
sex.config['PARAMETERS_LIST'].append('PETRO_RADIUS')
sex.config['PARAMETERS_LIST'].append('ISOAREA_IMAGE')
sex.config['PARAMETERS_LIST'].append('ALPHA_J2000')
sex.config['PARAMETERS_LIST'].append('DELTA_J2000')
sex.config['PARAMETERS_LIST'].append('FWHM_WORLD')
sex.config['PARAMETERS_LIST'].append('CLASS_STAR')
sex.config['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
sex.run(fname)
segmap = fits.open('check.fits')[0].data

df_cat = pd.read_table('py-sextractor.cat', delim_whitespace=True, header=16)
df_cat.columns = ['num','flux_best','fluxerr_best', 'x','y','flags',
               'fwhm_image', 'flux_iso','mag_isocor','mag_auto',
              'petro_radius','ISOAREA_IMAGE','ra','dec',
              'fwhm_world','class_star']

#selecao dos objetos que devem ser galaxias
df_cat = df_cat.ix[(df_cat['fwhm_image'] > 4.5) & (df_cat['mag_auto'] < -7)]
'''
================================================================================
Lendo as imagens, em todas as bandas, e gerando um dataframe para cada galaxia
utilizando astropy
Calculando o ceu em todas as bandas

ATUALIZAR NOME DA BANDA DE SEGMENTACAO
================================================================================
'''

df = pd.DataFrame()
df_sky = pd.DataFrame()

for i_object in range(len(df_cat)):
    window_size = 250
    filter_seg = 'rSDSS'
    ra = df_cat['ra']
    dec = df_cat['dec']
    image_r = fits.open('data/frame-r-002507-4-0226.fits')
    wcsys = wcs.WCS(header=image_r[0].header)
    y, x = wcsys.wcs_world2pix(ra, dec, 1)
    interval = (int(round(x[i_object] - window_size / 2)), int(round(x[i_object] + window_size / 2)),
                int(round(y[i_object] - window_size / 2)), int(round(y[i_object] + window_size / 2)))
    df = pd.DataFrame()
    df_sky = pd.DataFrame()
    seg_sex = segmap[interval[0]:interval[1], interval[2]:interval[3]]
    print('objeto: %d' %i_object)

    for i_gal in range(len(df_fit)):
        print('banda: %s' %df_fit['filter'][i_gal])
        f_sdss = fits.open('data/frame-%s-%s' %(df_fit['filter'][i_gal],
                                                 df_fit['name'][i_gal]))
        img = get_image(f_sdss)
#        a = int(interval[0])
        img_cut = img[interval[0]:interval[1], interval[2]:interval[3]]
        plt.figure(1)
        plt.clf()
        plt.imshow(100*np.log10(img_cut/255), cmap='spectral')
        plt.colorbar()
        band=df_fit['filter'][i_gal]
        nrows, ncols = img_cut.shape
        xx, yy = np.meshgrid( *np.ogrid[:ncols, :nrows] )
        table = np.column_stack(( xx.flatten(), yy.flatten(), img_cut.flatten() ))
        temp = pd.DataFrame(table, columns=['x','y',band])
        df = pd.concat([df,temp], axis=1)

        sky_r = fits.open('data/frame-%s-%s' %(df_fit['filter'][i_gal],
                                                 df_fit['name'][i_gal]))
        sky = get_image(sky_r)
        wcsys = wcs.WCS(header=sky_r[0].header)
        yc, xc = wcsys.wcs_world2pix(351.101, 14.737636, 1)
        delta_x = 85
        delta_y = 85
        interval_sky = (int(round(xc - delta_x / 2)), int(round(xc + delta_x / 2)), int(round(yc - delta_y / 2)),
                    int(round(yc + delta_y / 2)))
        img_sky = sky[interval_sky[0]:interval_sky[1], interval_sky[2]:interval_sky[3]]
        sky_nrows, sky_ncols = img_sky.shape
        xxc, yyc = np.meshgrid( *np.ogrid[:sky_ncols, :sky_nrows] )
        table_sky = np.column_stack(( xxc.flatten(), yyc.flatten(), img_sky.flatten() ))
        temp_sky = pd.DataFrame(table_sky, columns=['x','y',band])
        df_sky = pd.concat([df_sky,temp_sky], axis=1)

fim = time.time()
time_proc = fim - ini
print('')
print(bcolors.HEADER + 'tempo de processamento: %fs' %time_proc + bcolors.ENDC)
