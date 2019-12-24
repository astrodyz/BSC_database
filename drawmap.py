import os 
#os.environ['PROJ_LIB'] = '/Users/yuzhouyi/anaconda3/pkgs/proj4-5.2.0-h0a44026_1/share/proj'
# This line is special!!
import pandas as pd
import matplotlib.pyplot as plt
import datetime
from mpl_toolkits.basemap import Basemap
from matplotlib import path
import numpy as np
import sys

## Usage: python drawmap.py ra dec choice

def bluemarble_daynight1(date,lon, lat, scale):

    # Define Bluemarble and Nightshade objects
    fig, axes = plt.subplots(1, figsize=(16,16))
    m  = Basemap(projection='ortho', resolution= None, lat_0=lat[0], lon_0=lon[0], 
                area_thresh=None, ax=axes)
    bm = m.bluemarble(scale=scale)
    ns = m.nightshade(date, alpha=0.5)

    bm_rgb = bm.get_array()
    bm_ext = bm.get_extent()

    axes.cla()

    # Get the x and y index spacing
    x = np.linspace(bm_ext[0], bm_ext[1], bm_rgb.shape[1])
    y = np.linspace(bm_ext[2], bm_ext[3], bm_rgb.shape[0])

    # Define coordinates of the Bluemarble image
    x3d,y3d = np.meshgrid(x,y)
    pts     = np.hstack((x3d.flatten()[:,np.newaxis],y3d.flatten()[:,np.newaxis]))

    # Find which coordinates fall in Nightshade 
    # The following could be tidied up as there should only ever one polygon. Although
    # the length of ns.collections is 3? I'm sure there's a better way to do this.
    paths, polygons = [], []
    for i, polygons in enumerate(ns.collections):
        for j, paths in enumerate(polygons.get_paths()):
            #print j, i
            msk = paths.contains_points(pts)

    # Redefine mask
    msk        = np.reshape(msk,bm_rgb[:,:,0].shape)
    msk_s      = np.zeros(msk.shape)
    msk_s[~msk] = 1.

    # Smooth interface between Night and Day
    for s in range(int(bm_rgb.shape[1]/50)): # Make smoothing between day and night a function of Bluemarble resolution
        msk_s = 0.25 * (  np.vstack( (msk_s[-1,:            ], msk_s[:-1, :            ]) )  \
                        + np.vstack( (msk_s[1:,:            ], msk_s[0  , :            ]) )  \
                        + np.hstack( (msk_s[: ,0, np.newaxis], msk_s[:  , :-1          ]) )  \
                        + np.hstack( (msk_s[: ,1:           ], msk_s[:  , -1,np.newaxis]) ) )

    # Define new RGBA array
    bm_rgba = np.dstack((bm_rgb[:,:,0:3], msk_s))
    # Plot up Bluemarble Nightshade
    m    = Basemap(projection='ortho', resolution= None,  lat_0=lat[0], lon_0=lon[0], 
                   area_thresh=None, ax=axes)
    bm_n = m.warpimage('./earth_lights_lrg.jpg',scale=scale)
#from https://eoimages.gsfc.nasa.gov/images/imagerecords/55000/55167/earth_lights_lrg.jpg
    bm_d = m.imshow(bm_rgba)
    x, y = m(lon[0], lat[0])
    plt.plot(x, y, 'or', markersize=10)
    plt.title('Day/Night Map for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"), fontsize = 20)
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    plt.savefig('output/fig/0/f1.jpg')

    size = np.size(lon)
    xy = np.zeros([size,2])
    for i in range(size):
        x, y = m(lon[i], lat[i])
        xy[i,0] = x
        xy[i,1] = y
    print(xy)
    np.savetxt('output/fig/0/positionxy.dat', xy)


def bluemarble_daynight2(date,lon, lat, scale):

    # Define Bluemarble and Nightshade objects
    fig, axes = plt.subplots(1, figsize=(16,12))
    m  = Basemap(projection='cyl', resolution= None,  
                area_thresh=None, ax=axes)
    bm = m.bluemarble(scale=scale)
    ns = m.nightshade(date, alpha=0.5)

    bm_rgb = bm.get_array()
    bm_ext = bm.get_extent()

    axes.cla()

    # Get the x and y index spacing
    x = np.linspace(bm_ext[0], bm_ext[1], bm_rgb.shape[1])
    y = np.linspace(bm_ext[2], bm_ext[3], bm_rgb.shape[0])

    # Define coordinates of the Bluemarble image
    x3d,y3d = np.meshgrid(x,y)
    pts     = np.hstack((x3d.flatten()[:,np.newaxis],y3d.flatten()[:,np.newaxis]))

    # Find which coordinates fall in Nightshade 
    # The following could be tidied up as there should only ever one polygon. Although
    # the length of ns.collections is 3? I'm sure there's a better way to do this.
    paths, polygons = [], []
    for i, polygons in enumerate(ns.collections):
        for j, paths in enumerate(polygons.get_paths()):
            #print j, i
            msk = paths.contains_points(pts)

    # Redefine mask
    msk        = np.reshape(msk,bm_rgb[:,:,0].shape)
    msk_s      = np.zeros(msk.shape)
    msk_s[~msk] = 1.

    # Smooth interface between Night and Day
    for s in range(int(bm_rgb.shape[1]/50)): # Make smoothing between day and night a function of Bluemarble resolution
        msk_s = 0.25 * (  np.vstack( (msk_s[-1,:            ], msk_s[:-1, :            ]) )  \
                        + np.vstack( (msk_s[1:,:            ], msk_s[0  , :            ]) )  \
                        + np.hstack( (msk_s[: ,0, np.newaxis], msk_s[:  , :-1          ]) )  \
                        + np.hstack( (msk_s[: ,1:           ], msk_s[:  , -1,np.newaxis]) ) )

    # Define new RGBA array
    bm_rgba = np.dstack((bm_rgb[:,:,0:3], msk_s))
    # Plot up Bluemarble Nightshade
    m    = Basemap(projection='cyl', resolution= None, 
                   area_thresh=None, ax=axes)
    bm_n = m.warpimage('./earth_lights_lrg.jpg',scale=scale)
#from https://eoimages.gsfc.nasa.gov/images/imagerecords/55000/55167/earth_lights_lrg.jpg
    bm_d = m.imshow(bm_rgba)
    x, y = m(lon, lat)
    plt.plot(x, y, 'or', markersize=10)
    plt.title('Day/Night Map for %s (UTC)' % date.strftime("%d %b %Y %H:%M:%S"),fontsize = 15)
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    plt.savefig('f2.jpg')

#date = datetime.datetime.utcnow()
## date = datetime.datetime(2000,1,1,0,0,0) year month day hour minute second (UTC time)

info = pd.read_csv('input/input_time_position.csv')
#info = np.loadtxt('input/input_time_position.csv',delimiter = ',',skiprows = 1, 
#        usecols = np.r_[range(0, 6), range(21, 23)])



lon = info['longitude']
lat = info['latitude']
#print(lon,lat)
date = datetime.datetime(int(info['year'][0]),int(info['month'][0]),int(info['day'][0]),int(info['hour'][0]),int(info['minute'][0]),int(info['second'][0]))
#print(date)

choice = 0

if choice == 0:
    bluemarble_daynight1(date,lon, lat, 1.0)
else:
    bluemarble_daynight2(date,120, 30, 0.5)