#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 16:35:36 2019

@author: dyz
"""

import numpy as np 
import matplotlib.pyplot as plt
import os
import pandas as pd
import csv
import astropy
from subprocess import getoutput
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from astropy.time import Time
from datetime import datetime
from scipy.optimize import curve_fit 
import astropy.constants as c
from PyAstronomy import pyasl
from scipy import integrate 
import gc

def ra_convert(rah,ram,ras):
    return rah*15e0+ram*15e0/60e0+ras*15e0/3600

def de_convert(design,ded,dem,des):
    if design==0:
        return -(ded+dem/60e0+des/3600e0)
    elif design==1:
        return ded+dem/60e0+des/3600e0
    
def search_target_in_BSC_circle(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius):
    t=Time(datetime(year, month, day, hour, minute, second))
    tjd=(t.jd-2451544.5)/365
    BSC=pd.read_csv('data/BSC_catalog.csv',sep=',')
    target=[]
    c_ra=ra_convert(rah,ram,ras)
    c_de=de_convert(design,ded,dem,des)
    for i in range(0,len(BSC)):
        pmra=float(BSC['pmRA'][i])
        pmde=float(BSC['pmDE'][i])
        RAJ2000=float(BSC['RAJ2000'][i])
        DEJ2000=float(BSC['DEJ2000'][i])
        RA=RAJ2000+pmra*tjd
        DE=DEJ2000+pmde*tjd        
        if (RA-c_ra)**2+(DE-c_de)**2<=radius**2:
            target.append(i)
    return target

def search_target_in_BSC_square(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius):
    t=Time(datetime(year, month, day, hour, minute, second))
    tjd=(t.jd-2451544.5)/365
    BSC=pd.read_csv('data/BSC_catalog.csv',sep=',')
    target=[]
    c_ra=ra_convert(rah,ram,ras)
    c_de=de_convert(design,ded,dem,des)
    for i in range(0,len(BSC)):
        pmra=float(BSC['pmRA'][i])
        pmde=float(BSC['pmDE'][i])
        RAJ2000=float(BSC['RAJ2000'][i])
        DEJ2000=float(BSC['DEJ2000'][i])
        RA=RAJ2000+pmra*tjd
        DE=DEJ2000+pmde*tjd        
        if abs(RA-c_ra)<=radius/2 and abs(DE-c_de)<=radius/2:
            target.append(i)
    return target

def get_month(m):
    if m=='Jan':
        return 1
    elif m=='Feb':
        return 2
    elif m=='Mar':
        return 3
    elif m=='Apr':
        return 4
    elif m=='May':
        return 5
    elif m=='Jun':
        return 6
    elif m=='Jul':
        return 7
    elif m=='Aug':
        return 8
    elif m=='Sep':
        return 9
    elif m=='Oct':
        return 10
    elif m=='Nov':
        return 11
    elif m=='Dec':
        return 12
    
def convert_ra(ra):
    h=float(ra[0:2])
    m=float(ra[3:5])
    s=float(ra[6:11])
    return h*15+m*15./60+s*15/3600



def convert_dec(dec):
    sign=dec[0]
    d=float(dec[1:3])
    m=float(dec[4:6])
    s=float(dec[7:11])
    if sign=='-':
        return -(d*1+m*1./60+s*1/3600)
    else:
        return d*+m*1./60+s*1/3600

        

def search_target_in_Solarsys_circle(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius,flag2):
    t=Time(datetime(year, month, day, hour, minute, second))
    tjd=t.jd
    c_ra=ra_convert(rah,ram,ras)
    c_de=de_convert(design,ded,dem,des)    
    namelist=['Callisto',
              'Europa',
              'Ganymede',
              'Io',
              'Jupiter',
              'Mars',
              'Mercury',
              'Moon',
              'Neptune',
              'Saturn',
              'Sun',
              'Uranus',
              'Venus']
    
    ra=[]
    dec=[]
    vmag=[]
    name=[]
    for j in range(0,len(namelist)):
        table=pd.read_csv('Solar_system/'+namelist[j]+'/'+str(year)+'.csv',sep=',')
#       print(table['time'][0])
        tlist=[]
        if len(table)>0:
            for i in range(0,len(table)-1):
                y1=int(table['time'][i][1:5])
                m1=int(get_month(table['time'][i][6:9]))
                d1=int(table['time'][i][10:12])
                h1=int(table['time'][i][13:15])
                t1=Time(datetime(y1, m1, d1, h1, 0, 0))
                t1jd=t1.jd
                tlist.append(t1jd)
        

            for i in range(0,len(tlist)-1):
                if tjd>tlist[i] and tjd<=tlist[i+1]:
                    ra1=convert_ra(table['ra'][i])
                    ra2=convert_ra(table['ra'][i+1])
                    dec1=convert_dec(table['dec'][i])
                    dec2=convert_dec(table['dec'][i+1])
                    mag1=table['vmag'][i]
                    mag2=table['vmag'][i+1]
                    RA=(ra2-ra1)/(tlist[i+1]-tlist[i])*(tjd-tlist[i])+ra1
                    DE=(dec2-dec1)/(tlist[i+1]-tlist[i])*(tjd-tlist[i])+dec1
                    Vmag=(mag2-mag1)/(tlist[i+1]-tlist[i])*(tjd-tlist[i])+mag1
                    if (RA-c_ra)**2+(DE-c_de)**2<radius**2:
                        if Vmag<=flag2:                
                            ra.append(RA)
                            dec.append(DE)
                            vmag.append(Vmag)
                            name.append(namelist[j])
        
    return name,ra,dec,vmag
    
def search_target_in_Solarsys_square(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius,flag2):
    t=Time(datetime(year, month, day, hour, minute, second))
    tjd=t.jd
    c_ra=ra_convert(rah,ram,ras)
    c_de=de_convert(design,ded,dem,des)    
    namelist=['Callisto',
              'Europa',
              'Ganymede',
              'Io',
              'Jupiter',
              'Mars',
              'Mercury',
              'Moon',
              'Neptune',
              'Saturn',
              'Sun',
              'Uranus',
              'Venus']
    
    ra=[]
    dec=[]
    vmag=[]
    name=[]
    for j in range(0,len(namelist)):
        table=pd.read_csv('Solar_system/'+namelist[j]+'/'+str(year)+'.csv',sep=',')
#       print(table['time'][0])
        tlist=[]
        if len(table)>0:
            for i in range(0,len(table)-1):
                y1=int(table['time'][i][1:5])
                m1=int(get_month(table['time'][i][6:9]))
                d1=int(table['time'][i][10:12])
                h1=int(table['time'][i][13:15])
                t1=Time(datetime(y1, m1, d1, h1, 0, 0))
                t1jd=t1.jd
                tlist.append(t1jd)
        

            for i in range(0,len(tlist)-1):
                if tjd>tlist[i] and tjd<=tlist[i+1]:
                    ra1=convert_ra(table['ra'][i])
                    ra2=convert_ra(table['ra'][i+1])
                    dec1=convert_dec(table['dec'][i])
                    dec2=convert_dec(table['dec'][i+1])
                    mag1=table['vmag'][i]
                    mag2=table['vmag'][i+1]
                    RA=(ra2-ra1)/(tlist[i+1]-tlist[i])*(tjd-tlist[i])+ra1
                    DE=(dec2-dec1)/(tlist[i+1]-tlist[i])*(tjd-tlist[i])+dec1
                    Vmag=(mag2-mag1)/(tlist[i+1]-tlist[i])*(tjd-tlist[i])+mag1
                    if abs(RA-c_ra)<=radius/2 and abs(DE-c_de)<=radius/2:
                        if Vmag<=flag2:                
                            ra.append(RA)
                            dec.append(DE)
                            vmag.append(Vmag)
                            name.append(namelist[j])
        
    return name,ra,dec,vmag                    

    
def variable_search(i):
    v_table=pd.read_csv('result/variable/'+str(i)+'/variable.csv',sep=',')
    vsx_flag=v_table['VSX_flag'][0]
    v_flag=v_table['V_flag'][0]
    passband=v_table['Passband'][0]
    if passband[0]=='b':
        passband=passband[2:len(passband)]        
    max_mag=v_table['max'][0]
    min_mag=v_table['min'][0]
    delta_mag=v_table['delta'][0]
    if vsx_flag==0:
        return 0,'Not in VSX catalog',passband,max_mag,min_mag,delta_mag
    elif vsx_flag==1:
        if v_flag==0:
            return 1,'Variable',passband,max_mag,min_mag,delta_mag
        if v_flag==1:
            return 2,'Suspected variable',passband,max_mag,min_mag,delta_mag
        if v_flag==2:
            return 3,'Constant or non-existing',passband,max_mag,min_mag,delta_mag
        if v_flag==3:
            return 4,'Possible duplicate',passband,max_mag,min_mag,delta_mag   


        

                
def powerlaw(x,a,b):
    return a*x**b

'''
with open('input/input_time_position.csv','w') as f:
    s=['year',
       'month',
       'day',
       'hour',
       'minute',
       'second',
       'rah',
       'ram',
       'ras',
       'design',
       'ded',
       'dem',
       'des',
       'radius',
       'flag0',
       'flag1',
       'flag2',
       'flag3',
       'flag4',
       'low',
       'high',
       'longitude',
       'latitude',]
    
    writer=csv.DictWriter(f,fieldnames=s)
    writer.writeheader()
    writer.writerow({'year':2020,
                     'month':1,
                     'day':1,
                     'hour':0,
                     'minute':0,
                     'second':0,
                     'rah':3,
                     'ram':0,
                     'ras':0,
                     'design':0,
                     'ded':10,
                     'dem':0,
                     'des':0,
                     'radius':5,
                     'flag0':0,
                     'flag1':0,
                     'flag2':6.5,
                     'flag3':0,
                     'flag4':0,
                     'low':5000,
                     'high':6000,
                     'longitude':0,
                     'latitude':0,})
'''
  

#loat input file()
#BSClist=[]
input_f=pd.read_csv('input/input_time_position.csv',sep=',')
BSC=pd.read_csv('data/BSC_catalog.csv',sep=',')
VSX=pd.read_csv('data/VSX_catalog.csv',sep=',')
os.system('rm -rf output/fig')
os.system('rm -rf output/table')
for i in range(0,len(input_f)):
#flag0=0 BSC+solar system, flag0=1 BSC  
    os.system('mkdir output/fig')
    os.system('mkdir output/fig/'+str(i))
    os.system('mkdir output/table')
#    os.system('mkdir output/table/'+str(i))
    year=input_f['year'][i]
    month=input_f['month'][i]
    day=input_f['day'][i]
    hour=input_f['hour'][i]    
    minute=input_f['minute'][i]
    second=input_f['second'][i]
    rah=input_f['rah'][i]
    ram=input_f['ram'][i]
    ras=input_f['ras'][i]
    design=input_f['design'][i]
    ded=input_f['ded'][i]
    dem=input_f['dem'][i]
    des=input_f['des'][i]
    radius=input_f['radius'][i]
    flag0=input_f['flag0'][i]
    flag1=input_f['flag1'][i]
    flag2=input_f['flag2'][i]
    flag3=input_f['flag3'][i]
    flag4=input_f['flag4'][i]
    low=input_f['low'][i]
    high=input_f['high'][i]
    longitude=input_f['longitude'][i]
    latitude=input_f['latitude'][i]
    
        
    if flag0==0: 
        
        if flag4==0:
            BSClist=search_target_in_BSC_square(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius)
            Solarsyslist=search_target_in_Solarsys_square(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius,flag2)
        elif flag4==1:
            BSClist=search_target_in_BSC_circle(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius)
            Solarsyslist=search_target_in_Solarsys_circle(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius,flag2)
        
        with open('output/table/output_'+str(i)+'.csv','w') as f:
            s=['star_id',
               'name',               
               'RAJ2000',
               'DEJ2000',
               'Vmag_BSC',
#               'Bmag_BSC',
#               'Umag_BSC',
#               'Vmag_TIC',
#               'Jmag_TIC',
#               'Hmag_TIC',               
               'Variable_flag',
               'Passband',
               'max_mag',
               'min_mag',
               'delta_mag',
               's_flag',
               'e_flag',
               't_flag',
               'teff_flag',
               'teff',
               'flux',]
            writer=csv.DictWriter(f,fieldnames=s)
            writer.writeheader()
            #flag1=0 all, flag1=1 north, flag1=2 south
            if flag1==0:                         
                for j in range(0,len(BSClist)):   
                    os.system('mkdir output/fig/'+str(i)+'/'+str(BSClist[j]))
                    #flag2 vmag upper limit       
                    if BSC['Vmag'][BSClist[j]]<flag2:
                        
                        variable_table=pd.read_csv('result/'+str(BSClist[j])+'/variable_rewrite.csv')
                        spec_flag_table=pd.read_csv('result/'+str(BSClist[j])+'/spectrum_flag.csv')
                        
                        if flag3==0: #include variavle
                            #draw picture and calculate flux
                            if spec_flag_table['teff_flag'][0]==1:
                                spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                teff=spec_fit_table['temperature'][0] 
                                def blackbody(lam): 
                                    cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                    cof_teff=cof['temperature'][0]
                                    cof_cof=cof['coefficient'][0]

                                    from scipy.constants import h,k,c
                                    lam = 1e-10 * lam # convert to metres
                                    T=cof_teff
                                    return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                xx=np.linspace(3000,10000,701)
                                yy2=blackbody(xx)
                                fig=plt.figure(figsize=(7,6))
                                plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                flux=integrate.quad(blackbody,low,high)
                                obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                obv_phi=[]
                                for k in range(0,len(obv)):
                                    obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                plt.xlim(3000,10000)
                                plt.legend(loc ='best',fontsize=14)
                                plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        

                                plt.close()
                                
                                del xx,yy2
                                gc.collect()
                                os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            else:
                                flux=np.nan
                                teff=np.nan
                        
                        
                            if spec_flag_table['s_flag'][0]==1:
                                os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            if spec_flag_table['e_flag'][0]==1:
                                os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                            writer.writerow({'star_id':BSClist[j],
                                             'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                             'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                             'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                             'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                             'Variable_flag':variable_table['V_flag'][0],
                                             'Passband':variable_table['Passband'][0],
                                             'max_mag':variable_table['max_mag'][0],
                                             'min_mag':variable_table['min_mag'][0],
                                             'delta_mag':variable_table['delta_mag'][0],
                                             's_flag':spec_flag_table['s_flag'][0],
                                             'e_flag':spec_flag_table['e_flag'][0],
                                             't_flag':spec_flag_table['t_flag'][0],
                                             'teff_flag':spec_flag_table['teff_flag'][0],
                                             'teff':teff,
                                             'flux':flux[0],
                                             })
    
    
                        elif flag3==1: #exclude variable
                            if variable_table['V_flag'][0]==0:
                                if spec_flag_table['teff_flag'][0]==1:
                                    spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                    teff=spec_fit_table['temperature'][0] 
                                    def blackbody(lam): 
                                        cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                        cof_teff=cof['temperature'][0]
                                        cof_cof=cof['coefficient'][0]

                                        from scipy.constants import h,k,c
                                        lam = 1e-10 * lam # convert to metres
                                        T=cof_teff
                                        return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                    xx=np.linspace(3000,10000,701)
                                    yy2=blackbody(xx)
                                    fig=plt.figure(figsize=(7,6))
                                    plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                    flux=integrate.quad(blackbody,low,high)
                                    obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                    obv_phi=[]
                                    for k in range(0,len(obv)):
                                        obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                    plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                    plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                    plt.xlim(3000,10000)
                                    plt.legend(loc ='best',fontsize=14)
                                    plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                    plt.close()
                                    del xx,yy2
                                    gc.collect()
                                    os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                else:
                                    flux=np.nan
                                    teff=np.nan
                        
                        
                                if spec_flag_table['s_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                if spec_flag_table['e_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                writer.writerow({'star_id':BSClist[j],
                                                 'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                 'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                 'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                 'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                 'Variable_flag':variable_table['V_flag'][0],
                                                 'Passband':variable_table['Passband'][0],
                                                 'max_mag':variable_table['max_mag'][0],
                                                 'min_mag':variable_table['min_mag'][0],
                                                 'delta_mag':variable_table['delta_mag'][0],
                                                 's_flag':spec_flag_table['s_flag'][0],
                                                 'e_flag':spec_flag_table['e_flag'][0],
                                                 't_flag':spec_flag_table['t_flag'][0],
                                                 'teff_flag':spec_flag_table['teff_flag'][0],
                                                 'teff':teff,
                                                 'flux':flux[0],
                                                 })    

                            
                if len(Solarsyslist[0])>0:
                    for j in range(0,len(Solarsyslist[0])): #[0] name; [1] ra; [2] dec;[3] vmag
                        writer.writerow({'star_id':np.nan,
                                         'name':Solarsyslist[0][j],
                                         'RAJ2000':Solarsyslist[1][j],
                                         'DEJ2000':Solarsyslist[2][j],
                                         'Vmag_BSC':Solarsyslist[3][j],
                                         'Variable_flag':np.nan,
                                         'Passband':np.nan,
                                         'max_mag':np.nan,
                                         'min_mag':np.nan,
                                         'delta_mag':np.nan,
                                         's_flag':np.nan,
                                         'e_flag':np.nan,
                                         't_flag':np.nan,
                                         'teff_flag':np.nan,
                                         'teff':np.nan,
                                         'flux':np.nan,
                                         })    
                    
                                                        
    
            if flag1==1:
                for j in range(0,len(BSClist)):
                    if BSC['DEJ2000'][BSClist[j]]>0:
                        os.system('mkdir output/fig/'+str(i)+'/'+str(BSClist[j]))
                        #flag2 vmag upper limit       
                        if BSC['Vmag'][BSClist[j]]<flag2:
                            variable_table=pd.read_csv('result/'+str(BSClist[j])+'/variable_rewrite.csv')
                            spec_flag_table=pd.read_csv('result/'+str(BSClist[j])+'/spectrum_flag.csv')
                        
                            if flag3==0: #include variavle
                                #draw picture and calculate flux
                                if spec_flag_table['teff_flag'][0]==1:
                                    spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                    teff=spec_fit_table['temperature'][0] 
                                    def blackbody(lam): 
                                        cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                        cof_teff=cof['temperature'][0]
                                        cof_cof=cof['coefficient'][0]

                                        from scipy.constants import h,k,c
                                        lam = 1e-10 * lam # convert to metres
                                        T=cof_teff
                                        return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                    xx=np.linspace(3000,10000,701)
                                    yy2=blackbody(xx)
                                    fig=plt.figure(figsize=(7,6))
                                    plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                    flux=integrate.quad(blackbody,low,high)
                                    obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                    obv_phi=[]
                                    for k in range(0,len(obv)):
                                        obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                    plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                    plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                    plt.xlim(3000,10000)
                                    plt.legend(loc ='best',fontsize=14)
                                    plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                    plt.close()
                                    
                                    del xx,yy2
                                    gc.collect()
                                    os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                else:
                                    flux=np.nan
                                    teff=np.nan
                        
                        
                                if spec_flag_table['s_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                if spec_flag_table['e_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                writer.writerow({'star_id':BSClist[j],
                                                 'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                 'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                 'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                 'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                 'Variable_flag':variable_table['V_flag'][0],
                                                 'Passband':variable_table['Passband'][0],
                                                 'max_mag':variable_table['max_mag'][0],
                                                 'min_mag':variable_table['min_mag'][0],
                                                 'delta_mag':variable_table['delta_mag'][0],
                                                 's_flag':spec_flag_table['s_flag'][0],
                                                 'e_flag':spec_flag_table['e_flag'][0],
                                                 't_flag':spec_flag_table['t_flag'][0],
                                                 'teff_flag':spec_flag_table['teff_flag'][0],
                                                 'teff':teff,
                                                 'flux':flux[0],
                                                 })
    
    
                            elif flag3==1: #exclude variable
                                if variable_table['V_flag'][0]==0:
                                    if spec_flag_table['teff_flag'][0]==1:
                                        spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                        teff=spec_fit_table['temperature'][0] 
                                        def blackbody(lam): 
                                            cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                            cof_teff=cof['temperature'][0]
                                            cof_cof=cof['coefficient'][0]

                                            from scipy.constants import h,k,c
                                            lam = 1e-10 * lam # convert to metres
                                            T=cof_teff
                                            return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                        xx=np.linspace(3000,10000,701)
                                        yy2=blackbody(xx)
                                        fig=plt.figure(figsize=(7,6))
                                        plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                        flux=integrate.quad(blackbody,low,high)
                                        obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                        obv_phi=[]
                                        for k in range(0,len(obv)):
                                            obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                        plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                        plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                        plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                        plt.xlim(3000,10000)
                                        plt.legend(loc ='best',fontsize=14)
                                        plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                        plt.close()
                                        del xx,yy2
                                        gc.collect()
                                        os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    else:
                                        flux=np.nan
                                        teff=np.nan
                        
                            
                                    if spec_flag_table['s_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    if spec_flag_table['e_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                    writer.writerow({'star_id':BSClist[j],
                                                     'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                     'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                     'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                     'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                     'Variable_flag':variable_table['V_flag'][0],
                                                     'Passband':variable_table['Passband'][0],
                                                     'max_mag':variable_table['max_mag'][0],
                                                     'min_mag':variable_table['min_mag'][0],
                                                     'delta_mag':variable_table['delta_mag'][0],
                                                     's_flag':spec_flag_table['s_flag'][0],
                                                     'e_flag':spec_flag_table['e_flag'][0],
                                                     't_flag':spec_flag_table['t_flag'][0],
                                                     'teff_flag':spec_flag_table['teff_flag'][0],
                                                     'teff':teff,
                                                     'flux':flux[0],
                                                     })    

                            
                    if len(Solarsyslist[0])>0:
                        for j in range(0,len(Solarsyslist[0])): #[0] name; [1] ra; [2] dec;[3] vmag
                            writer.writerow({'star_id':np.nan,
                                             'name':Solarsyslist[0][j],
                                             'RAJ2000':Solarsyslist[1][j],
                                             'DEJ2000':Solarsyslist[2][j],
                                             'Vmag_BSC':Solarsyslist[3][j],
                                             'Variable_flag':np.nan,
                                             'Passband':np.nan,
                                             'max_mag':np.nan,
                                             'min_mag':np.nan,
                                             'delta_mag':np.nan,
                                             's_flag':np.nan,
                                             'e_flag':np.nan,
                                             't_flag':np.nan,
                                             'teff_flag':np.nan,
                                             'teff':np.nan,
                                             'flux':np.nan,
                                             })    


            if flag1==2:
                for j in range(0,len(BSClist)):
                    if BSC['DEJ2000'][BSClist[j]]<0:
                        os.system('mkdir output/fig/'+str(i)+'/'+str(BSClist[j]))
                        #flag2 vmag upper limit       
                        if BSC['Vmag'][BSClist[j]]<flag2:
                            variable_table=pd.read_csv('result/'+str(BSClist[j])+'/variable_rewrite.csv')
                            spec_flag_table=pd.read_csv('result/'+str(BSClist[j])+'/spectrum_flag.csv')
                        
                            if flag3==0: #include variavle
                                #draw picture and calculate flux
                                if spec_flag_table['teff_flag'][0]==1:
                                    spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                    teff=spec_fit_table['temperature'][0] 
                                    def blackbody(lam): 
                                        cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                        cof_teff=cof['temperature'][0]
                                        cof_cof=cof['coefficient'][0]

                                        from scipy.constants import h,k,c
                                        lam = 1e-10 * lam # convert to metres
                                        T=cof_teff
                                        return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                    xx=np.linspace(3000,10000,701)
                                    yy2=blackbody(xx)
                                    fig=plt.figure(figsize=(7,6))
                                    plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                    flux=integrate.quad(blackbody,low,high)
                                    obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                    obv_phi=[]
                                    for k in range(0,len(obv)):
                                        obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                    plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                    plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                    plt.xlim(3000,10000)
                                    plt.legend(loc ='best',fontsize=14)
                                    plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                    plt.close()
                                    del xx,yy2
                                    gc.collect()
                                    os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                else:
                                    flux=np.nan
                                    teff=np.nan
                        
                        
                                if spec_flag_table['s_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                if spec_flag_table['e_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                writer.writerow({'star_id':BSClist[j],
                                                 'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                 'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                 'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                 'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                 'Variable_flag':variable_table['V_flag'][0],
                                                 'Passband':variable_table['Passband'][0],
                                                 'max_mag':variable_table['max_mag'][0],
                                                 'min_mag':variable_table['min_mag'][0],
                                                 'delta_mag':variable_table['delta_mag'][0],
                                                 's_flag':spec_flag_table['s_flag'][0],
                                                 'e_flag':spec_flag_table['e_flag'][0],
                                                 't_flag':spec_flag_table['t_flag'][0],
                                                 'teff_flag':spec_flag_table['teff_flag'][0],
                                                 'teff':teff,
                                                 'flux':flux[0],
                                                 })
    
    
                            elif flag3==1: #exclude variable
                                if variable_table['V_flag'][0]==0:
                                    if spec_flag_table['teff_flag'][0]==1:
                                        spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                        teff=spec_fit_table['temperature'][0] 
                                        def blackbody(lam): 
                                            cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                            cof_teff=cof['temperature'][0]
                                            cof_cof=cof['coefficient'][0]

                                            from scipy.constants import h,k,c
                                            lam = 1e-10 * lam # convert to metres
                                            T=cof_teff
                                            return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                        xx=np.linspace(3000,10000,701)
                                        yy2=blackbody(xx)
                                        fig=plt.figure(figsize=(7,6))
                                        plt.plot(xx,yy2,'b-',label='Blackbody fit')
                                        
                                        flux=integrate.quad(blackbody,low,high)
                                        obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                        obv_phi=[]
                                        for k in range(0,len(obv)):
                                            obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
                                            
                                        plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                        plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                        plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                        plt.xlim(3000,10000)
                                        plt.legend(loc ='best',fontsize=14)
                                        plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                        plt.close()
                                        del xx,yy2
                                        gc.collect()
                                        os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    else:
                                        flux=np.nan
                                        teff=np.nan
                        
                            
                                    if spec_flag_table['s_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    if spec_flag_table['e_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                    writer.writerow({'star_id':BSClist[j],
                                                     'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                     'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                     'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                     'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                     'Variable_flag':variable_table['V_flag'][0],
                                                     'Passband':variable_table['Passband'][0],
                                                     'max_mag':variable_table['max_mag'][0],
                                                     'min_mag':variable_table['min_mag'][0],
                                                     'delta_mag':variable_table['delta_mag'][0],
                                                     's_flag':spec_flag_table['s_flag'][0],
                                                     'e_flag':spec_flag_table['e_flag'][0],
                                                     't_flag':spec_flag_table['t_flag'][0],
                                                     'teff_flag':spec_flag_table['teff_flag'][0],
                                                     'teff':teff,
                                                     'flux':flux[0],
                                                     })    

                            
                    if len(Solarsyslist[0])>0:
                        for j in range(0,len(Solarsyslist[0])): #[0] name; [1] ra; [2] dec;[3] vmag
                            writer.writerow({'star_id':np.nan,
                                             'name':Solarsyslist[0][j],
                                             'RAJ2000':Solarsyslist[1][j],
                                             'DEJ2000':Solarsyslist[2][j],
                                             'Vmag_BSC':Solarsyslist[3][j],
                                             'Variable_flag':np.nan,
                                             'Passband':np.nan,
                                             'max_mag':np.nan,
                                             'min_mag':np.nan,
                                             'delta_mag':np.nan,
                                             's_flag':np.nan,
                                             'e_flag':np.nan,
                                             't_flag':np.nan,
                                             'teff_flag':np.nan,
                                             'teff':np.nan,
                                             'flux':np.nan,
                                             })    

    if flag0==1: 
        
        if flag4==0:
            BSClist=search_target_in_BSC_square(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius)
#            Solarsyslist=search_target_in_Solarsys_square(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius,flag2)
        elif flag4==1:
            BSClist=search_target_in_BSC_circle(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius)
#            Solarsyslist=search_target_in_Solarsys_circle(year,month,day,hour,minute,second,rah,ram,ras,design,ded,dem,des,radius,flag2)
        
        with open('output/table/output_'+str(i)+'.csv','w') as f:
            s=['star_id',
               'name',               
               'RAJ2000',
               'DEJ2000',
               'Vmag_BSC',
#               'Bmag_BSC',
#               'Umag_BSC',
#               'Vmag_TIC',
#               'Jmag_TIC',
#               'Hmag_TIC',               
               'Variable_flag',
               'Passband',
               'max_mag',
               'min_mag',
               'delta_mag',
               's_flag',
               'e_flag',
               't_flag',
               'teff_flag',
               'teff',
               'flux',]
            writer=csv.DictWriter(f,fieldnames=s)
            writer.writeheader()
            #flag1=0 all, flag1=1 north, flag1=2 south
            if flag1==0:                         
                for j in range(0,len(BSClist)):
                    os.system('mkdir output/fig/'+str(i)+'/'+str(BSClist[j]))
                    #flag2 vmag upper limit       
                    if BSC['Vmag'][BSClist[j]]<flag2:
                        variable_table=pd.read_csv('result/'+str(BSClist[j])+'/variable_rewrite.csv')
                        spec_flag_table=pd.read_csv('result/'+str(BSClist[j])+'/spectrum_flag.csv')
                        
                        if flag3==0: #include variavle
                            #draw picture and calculate flux
                            if spec_flag_table['teff_flag'][0]==1:
                                spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                teff=spec_fit_table['temperature'][0] 
                                def blackbody(lam): 
                                    cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                    cof_teff=cof['temperature'][0]
                                    cof_cof=cof['coefficient'][0]

                                    from scipy.constants import h,k,c
                                    lam = 1e-10 * lam # convert to metres
                                    T=cof_teff
                                    return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                xx=np.linspace(3000,10000,701)
                                yy2=blackbody(xx)
                                fig=plt.figure(figsize=(7,6))
                                plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                flux=integrate.quad(blackbody,low,high)
                                obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                obv_phi=[]
                                for k in range(0,len(obv)):
                                    obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                plt.xlim(3000,10000)
                                plt.legend(loc ='best',fontsize=14)
                                plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                plt.close()
                                del xx,yy2
                                gc.collect()
                                os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            else:
                                flux=np.nan
                                teff=np.nan
                        
                        
                            if spec_flag_table['s_flag'][0]==1:
                                os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            if spec_flag_table['e_flag'][0]==1:
                                os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                            writer.writerow({'star_id':BSClist[j],
                                             'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                             'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                             'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                             'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                             'Variable_flag':variable_table['V_flag'][0],
                                             'Passband':variable_table['Passband'][0],
                                             'max_mag':variable_table['max_mag'][0],
                                             'min_mag':variable_table['min_mag'][0],
                                             'delta_mag':variable_table['delta_mag'][0],
                                             's_flag':spec_flag_table['s_flag'][0],
                                             'e_flag':spec_flag_table['e_flag'][0],
                                             't_flag':spec_flag_table['t_flag'][0],
                                             'teff_flag':spec_flag_table['teff_flag'][0],
                                             'teff':teff,
                                             'flux':flux[0],
                                             })
    
    
                        elif flag3==1: #exclude variable
                            if variable_table['V_flag'][0]==0:
                                if spec_flag_table['teff_flag'][0]==1:
                                    spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                    teff=spec_fit_table['temperature'][0] 
                                    def blackbody(lam): 
                                        cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                        cof_teff=cof['temperature'][0]
                                        cof_cof=cof['coefficient'][0]

                                        from scipy.constants import h,k,c
                                        lam = 1e-10 * lam # convert to metres
                                        T=cof_teff
                                        return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                    xx=np.linspace(3000,10000,701)
                                    yy2=blackbody(xx)
                                    fig=plt.figure(figsize=(7,6))
                                    plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                    flux=integrate.quad(blackbody,low,high)
                                    obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                    obv_phi=[]
                                    for k in range(0,len(obv)):
                                        obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                    plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                    plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                    plt.xlim(3000,10000)
                                    plt.legend(loc ='best',fontsize=14)
                                    plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                    plt.close()
                                    del xx,yy2
                                    gc.collect()
                                    os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                else:
                                    flux=np.nan
                                    teff=np.nan
                        
                        
                                if spec_flag_table['s_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                if spec_flag_table['e_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                writer.writerow({'star_id':BSClist[j],
                                                 'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                 'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                 'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                 'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                 'Variable_flag':variable_table['V_flag'][0],
                                                 'Passband':variable_table['Passband'][0],
                                                 'max_mag':variable_table['max_mag'][0],
                                                 'min_mag':variable_table['min_mag'][0],
                                                 'delta_mag':variable_table['delta_mag'][0],
                                                 's_flag':spec_flag_table['s_flag'][0],
                                                 'e_flag':spec_flag_table['e_flag'][0],
                                                 't_flag':spec_flag_table['t_flag'][0],
                                                 'teff_flag':spec_flag_table['teff_flag'][0],
                                                 'teff':teff,
                                                 'flux':flux[0],
                                                 })    

                '''                            
                if len(Solarsyslist[0])>0:
                    for j in range(0,len(Solarsyslist[0])): #[0] name; [1] ra; [2] dec;[3] vmag
                        writer.writerow({'star_id':np.nan,
                                         'name':Solarsyslist[0][j],
                                         'RAJ2000':Solarsyslist[1][j],
                                         'DEJ2000':Solarsyslist[2][j],
                                         'Vmag_BSC':Solarsyslist[3][j],
                                         'Variable_flag':np.nan,
                                         'Passband':np.nan,
                                         'max_mag':np.nan,
                                         'min_mag':np.nan,
                                         'delta_mag':np.nan,
                                         's_flag':np.nan,
                                         'e_flag':np.nan,
                                         't_flag':np.nan,
                                         'teff_flag':np.nan,
                                         'teff':np.nan,
                                         'flux':np.nan,
                                         })    
                    '''                    
                                                        
    
            if flag1==1:
                for j in range(0,len(BSClist)):
                    if BSC['DEJ2000'][BSClist[j]]>0:
                        os.system('mkdir output/fig/'+str(i)+'/'+str(BSClist[j]))
                        #flag2 vmag upper limit       
                        if BSC['Vmag'][BSClist[j]]<flag2:
                            variable_table=pd.read_csv('result/'+str(BSClist[j])+'/variable_rewrite.csv')
                            spec_flag_table=pd.read_csv('result/'+str(BSClist[j])+'/spectrum_flag.csv')
                        
                            if flag3==0: #include variavle
                                #draw picture and calculate flux
                                if spec_flag_table['teff_flag'][0]==1:
                                    spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                    teff=spec_fit_table['temperature'][0] 
                                    def blackbody(lam): 
                                        cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                        cof_teff=cof['temperature'][0]
                                        cof_cof=cof['coefficient'][0]

                                        from scipy.constants import h,k,c
                                        lam = 1e-10 * lam # convert to metres
                                        T=cof_teff
                                        return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                    xx=np.linspace(3000,10000,701)
                                    yy2=blackbody(xx)
                                    fig=plt.figure(figsize=(7,6))
                                    plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                    flux=integrate.quad(blackbody,low,high)
                                    obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                    obv_phi=[]
                                    for k in range(0,len(obv)):
                                        obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                    plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                    plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                    plt.xlim(3000,10000)
                                    plt.legend(loc ='best',fontsize=14)
                                    plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                    plt.close()
                                    del xx,yy2
                                    gc.collect()
                                    os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                else:
                                    flux=np.nan
                                    teff=np.nan
                        
                        
                                if spec_flag_table['s_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                if spec_flag_table['e_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                writer.writerow({'star_id':BSClist[j],
                                                 'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                 'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                 'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                 'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                 'Variable_flag':variable_table['V_flag'][0],
                                                 'Passband':variable_table['Passband'][0],
                                                 'max_mag':variable_table['max_mag'][0],
                                                 'min_mag':variable_table['min_mag'][0],
                                                 'delta_mag':variable_table['delta_mag'][0],
                                                 's_flag':spec_flag_table['s_flag'][0],
                                                 'e_flag':spec_flag_table['e_flag'][0],
                                                 't_flag':spec_flag_table['t_flag'][0],
                                                 'teff_flag':spec_flag_table['teff_flag'][0],
                                                 'teff':teff,
                                                 'flux':flux[0],
                                                 })
    
    
                            elif flag3==1: #exclude variable
                                if variable_table['V_flag'][0]==0:
                                    if spec_flag_table['teff_flag'][0]==1:
                                        spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                        teff=spec_fit_table['temperature'][0] 
                                        def blackbody(lam): 
                                            cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                            cof_teff=cof['temperature'][0]
                                            cof_cof=cof['coefficient'][0]

                                            from scipy.constants import h,k,c
                                            lam = 1e-10 * lam # convert to metres
                                            T=cof_teff
                                            return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                        xx=np.linspace(3000,10000,701)
                                        yy2=blackbody(xx)
                                        fig=plt.figure(figsize=(7,6))
                                        plt.plot(xx,yy2,'b-',label='Blackbody fit')
                                        
                                        flux=integrate.quad(blackbody,low,high)
                                        obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                        obv_phi=[]
                                        for k in range(0,len(obv)):
                                            obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                        plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                        plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                        plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                        plt.xlim(3000,10000)
                                        plt.legend(loc ='best',fontsize=14)
                                        plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                        plt.close()
                                        del xx,yy2
                                        gc.collect()
                                        os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    else:
                                        flux=np.nan
                                        teff=np.nan
                        
                            
                                    if spec_flag_table['s_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    if spec_flag_table['e_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                    writer.writerow({'star_id':BSClist[j],
                                                     'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                     'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                     'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                     'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                     'Variable_flag':variable_table['V_flag'][0],
                                                     'Passband':variable_table['Passband'][0],
                                                     'max_mag':variable_table['max_mag'][0],
                                                     'min_mag':variable_table['min_mag'][0],
                                                     'delta_mag':variable_table['delta_mag'][0],
                                                     's_flag':spec_flag_table['s_flag'][0],
                                                     'e_flag':spec_flag_table['e_flag'][0],
                                                     't_flag':spec_flag_table['t_flag'][0],
                                                     'teff_flag':spec_flag_table['teff_flag'][0],
                                                     'teff':teff,
                                                     'flux':flux[0],
                                                     })    

                    '''        
                    if len(Solarsyslist[0])>0:
                        for j in range(0,len(Solarsyslist[0])): #[0] name; [1] ra; [2] dec;[3] vmag
                            writer.writerow({'star_id':np.nan,
                                             'name':Solarsyslist[0][j],
                                             'RAJ2000':Solarsyslist[1][j],
                                             'DEJ2000':Solarsyslist[2][j],
                                             'Vmag_BSC':Solarsyslist[3][j],
                                             'Variable_flag':np.nan,
                                             'Passband':np.nan,
                                             'max_mag':np.nan,
                                             'min_mag':np.nan,
                                             'delta_mag':np.nan,
                                             's_flag':np.nan,
                                             'e_flag':np.nan,
                                             't_flag':np.nan,
                                             'teff_flag':np.nan,
                                             'teff':np.nan,
                                             'flux':np.nan,
                                             })    
                    '''

            if flag1==2:
                for j in range(0,len(BSClist)):
                    if BSC['DEJ2000'][BSClist[j]]<0:
                        os.system('mkdir output/fig/'+str(i)+'/'+str(BSClist[j]))
                        #flag2 vmag upper limit       
                        if BSC['Vmag'][BSClist[j]]<flag2:
                            variable_table=pd.read_csv('result/'+str(BSClist[j])+'/variable_rewrite.csv')
                            spec_flag_table=pd.read_csv('result/'+str(BSClist[j])+'/spectrum_flag.csv')
                        
                            if flag3==0: #include variavle
                                #draw picture and calculate flux
                                if spec_flag_table['teff_flag'][0]==1:
                                    spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                    teff=spec_fit_table['temperature'][0] 
                                    def blackbody(lam): 
                                        cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                        cof_teff=cof['temperature'][0]
                                        cof_cof=cof['coefficient'][0]

                                        from scipy.constants import h,k,c
                                        lam = 1e-10 * lam # convert to metres
                                        T=cof_teff
                                        return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                    xx=np.linspace(3000,10000,701)
                                    yy2=blackbody(xx)
                                    fig=plt.figure(figsize=(7,6))
                                    plt.plot(xx,yy2,'b-',label='Blackbody fit')
                               
                                    flux=integrate.quad(blackbody,low,high)
                                    obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                    obv_phi=[]
                                    for k in range(0,len(obv)):
                                        obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
            
                                    plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                    plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                    plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                    plt.xlim(3000,10000)
                                    plt.legend(loc ='best',fontsize=14)
                                    plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                    plt.close()
                                    del xx,yy2
                                    gc.collect()
                                    os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                else:
                                    flux=np.nan
                                    teff=np.nan
                        
                        
                                if spec_flag_table['s_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                if spec_flag_table['e_flag'][0]==1:
                                    os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                writer.writerow({'star_id':BSClist[j],
                                                 'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                 'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                 'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                 'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                 'Variable_flag':variable_table['V_flag'][0],
                                                 'Passband':variable_table['Passband'][0],
                                                 'max_mag':variable_table['max_mag'][0],
                                                 'min_mag':variable_table['min_mag'][0],
                                                 'delta_mag':variable_table['delta_mag'][0],
                                                 's_flag':spec_flag_table['s_flag'][0],
                                                 'e_flag':spec_flag_table['e_flag'][0],
                                                 't_flag':spec_flag_table['t_flag'][0],
                                                 'teff_flag':spec_flag_table['teff_flag'][0],
                                                 'teff':teff,
                                                 'flux':flux[0],
                                                 })
    
    
                            elif flag3==1: #exclude variable
                                if variable_table['V_flag'][0]==0:
                                    if spec_flag_table['teff_flag'][0]==1:
                                        spec_fit_table=pd.read_csv('result/'+str(BSClist[j])+'/spectra_fit.csv')
                                        teff=spec_fit_table['temperature'][0] 
                                        def blackbody(lam): 
                                            cof=pd.read_csv('result/'+str(BSClist[j])+'/'+'spectra_fit.csv')
                                            cof_teff=cof['temperature'][0]
                                            cof_cof=cof['coefficient'][0]

                                            from scipy.constants import h,k,c
                                            lam = 1e-10 * lam # convert to metres
                                            T=cof_teff
                                            return cof_cof*2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))/(h*c/lam)*1e-7
                            
                                        xx=np.linspace(3000,10000,701)
                                        yy2=blackbody(xx)
                                        fig=plt.figure(figsize=(7,6))
                                        plt.plot(xx,yy2,'b-',label='Blackbody fit')
                                        
                                        flux=integrate.quad(blackbody,low,high)
                                        obv=pd.read_csv('result/'+str(BSClist[j])+'/'+'observation_data.csv')
                                        obv_phi=[]
                                        for k in range(0,len(obv)):
                                            obv_phi.append(obv['flux'][k]/(c.h.value*c.c.value/(obv['lamda'][k]*1e-10))*1e-7)
                                            
                                        plt.plot(obv['lamda'],obv_phi,'r.',label='Observation data')      
                                        plt.xlabel('Wavelength [$\AA$]',fontsize=14)
                                        plt.ylabel('Flux [photons cm$^{-1}$ s$^{-1}$ $\AA$ $^{-1}$]',fontsize=14)
                                        plt.xlim(3000,10000)
                                        plt.legend(loc ='best',fontsize=14)
                                        plt.savefig('output/fig/'+str(i)+'/'+str(BSClist[j])+'/spectra_fit.jpg',dpi=500)        
                                        plt.close()
                                        del xx,yy2
                                        gc.collect()
                                        os.system('cp result/'+str(BSClist[j])+'/'+'observation_data.csv '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    else:
                                        flux=np.nan
                                        teff=np.nan
                        
                            
                                    if spec_flag_table['s_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_sophie.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                                    if spec_flag_table['e_flag'][0]==1:
                                        os.system('cp result/'+str(BSClist[j])+'/'+'spectra_elodies.jpg '+'output/fig/'+str(i)+'/'+str(BSClist[j]))
                                
                            #flag3=0 include Variable, flag3=1 exclude Variable
                                    writer.writerow({'star_id':BSClist[j],
                                                     'name':'HD'+str(BSC['HD'][BSClist[j]]),
                                                     'RAJ2000':BSC['RAJ2000'][BSClist[j]],
                                                     'DEJ2000':BSC['DEJ2000'][BSClist[j]],
                                                     'Vmag_BSC':BSC['Vmag'][BSClist[j]],
                                                     'Variable_flag':variable_table['V_flag'][0],
                                                     'Passband':variable_table['Passband'][0],
                                                     'max_mag':variable_table['max_mag'][0],
                                                     'min_mag':variable_table['min_mag'][0],
                                                     'delta_mag':variable_table['delta_mag'][0],
                                                     's_flag':spec_flag_table['s_flag'][0],
                                                     'e_flag':spec_flag_table['e_flag'][0],
                                                     't_flag':spec_flag_table['t_flag'][0],
                                                     'teff_flag':spec_flag_table['teff_flag'][0],
                                                     'teff':teff,
                                                     'flux':flux[0],
                                                     })    

                            
                    '''        
                    if len(Solarsyslist[0])>0:
                        for j in range(0,len(Solarsyslist[0])): #[0] name; [1] ra; [2] dec;[3] vmag
                            writer.writerow({'star_id':np.nan,
                                             'name':Solarsyslist[0][j],
                                             'RAJ2000':Solarsyslist[1][j],
                                             'DEJ2000':Solarsyslist[2][j],
                                             'Vmag_BSC':Solarsyslist[3][j],
                                             'Variable_flag':np.nan,
                                             'Passband':np.nan,
                                             'max_mag':np.nan,
                                             'min_mag':np.nan,
                                             'delta_mag':np.nan,
                                             's_flag':np.nan,
                                             'e_flag':np.nan,
                                             't_flag':np.nan,
                                             'teff_flag':np.nan,
                                             'teff':np.nan,
                                             'flux':np.nan,
                                             })    
                    '''

    '''
def spectrum(i,low,high,l):   
    f1=open('result/spectrum/'+str(i)+'/'+'sophie_flag.dat','r')
    temp1=f1.readline()
    s_flag=temp1[0]
    f1.close()
    
    f2=open('result/spectrum/'+str(i)+'/'+'elodies_flag.dat','r')
    temp2=f2.readline()
    e_flag=temp2[0]
    f2.close()
    
    f3=open('result/spectrum/'+str(i)+'/'+'tic_flag.dat','r')
    temp3=f3.readline()
    t_flag=temp3[0]
    f3.close()
    
    zero_flux=999.5*1000
    if s_flag=='0':
        s_flux='Not in Sophie'
    elif s_flag=='1':
        print('Find it in Sophie')
        temp1=temp1.strip('\n')
        s_filename=temp1[2:len(temp1)+1]
        data1=np.loadtxt('result/spectrum/'+str(i)+'/'+s_filename,skiprows=4)
        wavelength=[]
        flux=[]
        v_r_flux=0
        r_flux=0
        for k in range(0,len(data1)):
            wavelength.append(data1[k][0])
            flux.append(data1[k][1])
                
        if 6000<=max(wavelength) and 5000>=min(wavelength):
            delta=wavelength[1]-wavelength[0]
            for k in range(0,len(wavelength)-1):
                if 6000>wavelength[k] and 6000<wavelength[k+1]:
                    high_id=k
                if 5000>wavelength[k] and 5000<wavelength[k+1]:
                    low_id=k
            for k in range(low_id,high_id):
                v_r_flux+=flux[k]*delta
        
        flux_ratio=zero_flux/v_r_flux*10**(-0.4*BSC['Vmag'][i])    

        
        for k in range(0,len(data1)):
            flux[k]=flux[k]*flux_ratio       
            
        bin_length=int(len(data1)/3)
        temp=0
        count=0
        for k in range(0,bin_length):
            if flux[k]==flux[k]:
                temp+=flux[k]
                count+=1
        y1=temp/count    
    
        temp=0
        count=0
        for k in range(bin_length,2*bin_length):
            if flux[k]==flux[k]:
                temp+=flux[k]
                count+=1
        y2=temp/count    

        temp=0
        count=0
        
        for k in range(2*bin_length,len(data1)):
            if flux[k]==flux[k]:
                temp+=flux[k]
                count+=1
        y3=temp/count    

        x1=np.mean(wavelength[0:bin_length])
        x2=np.mean(wavelength[bin_length:2*bin_length])
        x3=np.mean(wavelength[2*bin_length:len(data1)])
        
        x=[x1,x2,x3]
        y=[y1,y2,y3]
        popt, pcov = curve_fit(powerlaw, x, y)
        a=popt[0]
        b=popt[1]
        xx=np.linspace(4000,10000,10000)
        yy=powerlaw(xx,a,b)
        fig=plt.figure(figsize=(7,6))
        plt.loglog(wavelength,flux,'r-',label='Elodies')
        plt.loglog(xx,yy,'k--',label='power law fit')
        plt.xlabel('Wavelength [A]',fontsize=14)
        plt.ylabel('Flux [photons cm$^{-1}$ s^${-1}$ A$^{-1}$]',fontsize=14)
        plt.legend(loc ='best')
        plt.savefig('output/fig/'+str(l)+'/HD'+str(BSC['HD'][i])+'_spectrum_Sophie.eps')
        
        flux_print=0
        if high<=10000 and low>=4000:
            if high<=max(wavelength) and low>=min(wavelength):
                delta=wavelength[1]-wavelength[0]
                for k in range(0,len(wavelength)-1):
                    if high>wavelength[k] and high<wavelength[k+1]:
                        high_id=k
                    if low>wavelength[k] and low<wavelength[k+1]:
                        low_id=k
                        
                for k in range(low_id,high_id):
                    flux_print+=flux[k]*delta            
            
            else:
                flux_print=a/(b+1)*(high**(b+1)-low**(b+1))
                
            s_flux=str(flux_print)
        else:
            s_flux='out of range'
                    
        
    
    if e_flag=='0':
        e_flux='Not in Elodies'
    elif e_flag=='1':
        temp2=temp2.strip('\n')
        e_filename=temp2[2:len(temp2)+1]
        data2=np.loadtxt('result/spectrum/'+str(i)+'/'+e_filename,skiprows=4)
        wavelength=[]
        flux=[]
        v_r_flux=0
        r_flux=0
        for k in range(0,len(data2)):
            wavelength.append(data2[k][0])
            flux.append(data2[k][1])
                
        if 6000<=max(wavelength) and 5000>=min(wavelength):
            delta=wavelength[1]-wavelength[0]
            for k in range(0,len(wavelength)-1):
                if 6000>wavelength[k] and 6000<wavelength[k+1]:
                    high_id=k
                if 5000>wavelength[k] and 5000<wavelength[k+1]:
                    low_id=k
            for k in range(low_id,high_id):
                v_r_flux+=flux[k]*delta
        flux_ratio=zero_flux/v_r_flux*10**(-0.4*BSC['Vmag'][i])    


        for k in range(0,len(data2)):
            flux[k]=flux[k]*flux_ratio
        
        bin_length=int(len(data2)/3)
    
        temp=0
        count=0
        for k in range(0,bin_length):
            if flux[k]==flux[k]:
                temp+=flux[k]
                count+=1
        y1=temp/count    
    
        temp=0
        count=0
        for k in range(bin_length,2*bin_length):
            if flux[k]==flux[k]:
                temp+=flux[k]
                count+=1
    
        y2=temp/count    

        temp=0
        count=0
        for k in range(2*bin_length,len(data2)):
            if flux[k]==flux[k]:
                temp+=flux[k]
                count+=1
    
        y3=temp/count    
    
        x1=np.mean(wavelength[0:bin_length])
        x2=np.mean(wavelength[bin_length:2*bin_length])
        x3=np.mean(wavelength[2*bin_length:len(data2)])
    
        x=[x1,x2,x3]
        y=[y1,y2,y3]
        popt, pcov = curve_fit(powerlaw, x, y)
        a=popt[0]
        b=popt[1]
        xx=np.linspace(4000,10000,10000)
        yy=powerlaw(xx,a,b)
        fig=plt.figure(figsize=(7,6))
        plt.loglog(wavelength,flux,'r-',label='Elodies')
        plt.loglog(xx,yy,'k--',label='power law fit')
        plt.xlabel('Wavelength [A]',fontsize=14)
        plt.ylabel('Flux [photons cm$^{-1}$ s^${-1}$ A$^{-1}$]',fontsize=14)
        plt.legend(loc ='best')
        plt.savefig('output/fig/'+str(l)+'/HD'+str(BSC['HD'][i])+'_spectrum_Elodies.eps')
        
        flux_print=0
        if high<=10000 and low>=4000:
            if high<=max(wavelength) and low>=min(wavelength):
                delta=wavelength[1]-wavelength[0]
                for k in range(0,len(wavelength)-1):
                    if high>wavelength[k] and high<wavelength[k+1]:
                        high_id=k
                        
                    if low>wavelength[k] and low<wavelength[k+1]:
                        low_id=k
                        
                for k in range(low_id,high_id):
                    flux_print+=flux[k]*delta            
            
            else:
                flux_print=a/(b+1)*(high**(b+1)-low**(b+1))
                
            e_flux=str(flux_print)
            
        else:
            e_flux='out of range'


    
    if t_flag=='0':
        t_flux='Not in TIC'
    elif t_flag=='1':
        data3=pd.read_csv('result/spectrum/'+str(i)+'/tic.csv',sep=',')
        if data3['Teff'][0]!='--':
            temperature=float(data3['Teff'][0])*u.K
            w_value=np.linspace(low,high,1000)
            w=w_value*u.AA
            flux_lam=blackbody_lambda(w, temperature)
            f_value=flux_lam.value


        
            v_r_flux=0
            r_flux=0
            v_wavelength=np.linspace(5000,6000,1000)
            delta=v_wavelength[1]-v_wavelength[0]
            for k in range(0,len(v_wavelength)):
                flux_lam=blackbody_lambda(v_wavelength[k]*u.AA, temperature)
                w_temp=v_wavelength[k]
                f=c.c.value/w_temp*1e10
                e=f*c.h.value
                flux_temp=flux_lam.value*1e-7/e
                v_r_flux+=flux_temp*delta
                
            flux_ratio=zero_flux/v_r_flux*10**(-0.4*BSC['Vmag'][i])
    
            if high<=10000 and low>=4000: 
                v_wavelength=np.linspace(low,high,1000)
                delta=v_wavelength[1]-v_wavelength[0]
                for k in range(0,len(v_wavelength)):
                    flux_lam=blackbody_lambda(v_wavelength[k]*u.AA, temperature)
                    w_temp=v_wavelength[k]
                    f=c.c.value/w_temp*1e10
                    e=f*c.h.value
                    flux_temp=flux_lam.value*1e-7/e
                    r_flux+=flux_temp*delta
            
                flux=r_flux*flux_ratio
        
                t_flux=str(flux)
            else:
                t_flux='out of range'
    
            temperature=float(data3['Teff'][0])*u.K
            w_value=np.linspace(4000,10000,10000)
            delta=w_value[1]-w_value[0]
            w=w_value*u.AA
            flux_lam=blackbody_lambda(w, temperature)
            f_value=flux_lam.value
            fig=plt.figure(figsize=(7,6))
            plt.loglog(w_value,f_value,'r-',label='TIC')

            plt.xlabel('Wavelength [A]',fontsize=14)
            plt.ylabel('relative flux',fontsize=14)
            plt.legend(loc ='best')
            plt.savefig('output/fig/'+str(l)+'/HD'+str(BSC['HD'][i])+'_spectrum_TIC.eps')
        else:
            t_flux='No temperature'

    return s_flux,e_flux,t_flux
'''