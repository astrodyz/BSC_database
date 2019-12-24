from astropy.table import Table
import numpy as np
from photutils.datasets import (make_random_gaussians_table,
                                make_gaussian_sources_image,
                               make_noise_image)
import matplotlib.pyplot as plt
import os
import re

def drawpic(filename, data, ra, dec, radius, number):
    sigma_psf = 2.5
    sources = Table()
    size = np.size(data[:,0])
    x = ra
    y = dec
    x1 = x - radius
    x2 = x + radius
    y1 = y - radius
    y2 = y + radius
    sources['x_mean'] = (data[:,0] - x1) / radius / 2. * 256.
    sources['y_mean'] = (data[:,1] - y1) / radius / 2. * 256.
    sources['x_stddev'] = sigma_psf*np.ones(size)
    sources['y_stddev'] = sources['x_stddev']
    sources['theta'] = np.zeros(size)
    sources['flux'] = data[:,2] / np.min(data[:,2]) * 5000
    tshape = (256, 256)
    image = (make_gaussian_sources_image(tshape, sources) + \
            make_noise_image(tshape, distribution='poisson', mean=10.,
                         random_state=12) + \
            make_noise_image(tshape, distribution='gaussian', mean=0.,
                           stddev=10., random_state=12))

    plt.imshow(image, cmap='gray', extent = [x1, x2, y1, y2], interpolation='nearest',
           origin='lower') 
    plt.xlabel('RA(°)')
    plt.ylabel('DEC(°)')
    for i in range(size):
        plt.text(data[i,0],data[i,1],np.str(i), fontsize=12, color = 'w')
    plt.savefig('output/fig/' + np.str(number) + '/'+str(number)+'.jpg')
    

def sgn(i):
    if i == 0:
        return -1
    else:
        return 0



filename = os.popen('ls output/table/output_*.csv').read()
filename = filename.split('\n')
size = np.size(filename) - 1

info = np.loadtxt('input/input_time_position.csv',delimiter = ',',skiprows = 1, 
        usecols = np.r_[range(6, 14)])

if size==1:
    number = re.findall(r"\d+",filename[0])
    number = int(number[0])
    data = np.loadtxt(filename[0], delimiter = ',',skiprows = 1, usecols = [2,3,15])
    ra = (info[0] + info[1] / 60. + info[2] / 3600) * 15.0
    dec = sgn(info[3]) * (info[4] + info[5] / 60. + info[6] / 3600. )
    radius = info[7]

    drawpic(filename[0], data, ra, dec, radius, number)

else:    
    for i in range(size):
        number = re.findall(r"\d+",filename[i])
        number = int(number[0])
        data = np.loadtxt(filename[i], delimiter = ',',skiprows = 1, usecols = [2,3,15])
        ra = (info[number,0] + info[number,1] / 60. + info[number,2] / 3600) * 15.0
        dec = sgn(info[number,3]) * (info[number,4] + info[number,5] / 60. + info[number,6] / 3600. )
        radius = info[number,7]

        drawpic(filename[i], data, ra, dec, radius, number)





