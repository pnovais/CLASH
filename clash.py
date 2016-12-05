
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

df_fit = pd.read_table('data/arquivo_fits.txt', delim_whitespace=True)
df_fit.columns = ['filter','name']


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

# Lendo os nomes das imagens, em todas as bandas
# Para cada banda, iremos fazer o recorte e, posteriormente, uni-las em um unico
# dataframe pandas
#window_size = 150

fname = 'data/frame-r-002507-4-0226.fits'
sex = SExtractor()
sex.config['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
sex.run(fname)
segmap = fits.open('check.fits')[0].data

fim = time.time()
time_proc = fim - ini
print('')
print(bcolors.HEADER + 'tempo de processamento: %fs' %time_proc + bcolors.ENDC)
