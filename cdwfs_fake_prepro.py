#!/usr/bin/env python2
# -*- coding: utf-8 -*-
 
#%%  
import numpy as np
import os
import os.path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

from scipy.interpolate import interp1d
from scipy import integrate
from scipy.integrate import quad
from scipy.stats import chisquare
from scipy.optimize import curve_fit

from astropy.io import fits
from astropy import stats as sts
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.cosmology import WMAP9 as cosmo
from astropy.stats import sigma_clip

from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
 
from pylab import *

from shapely.geometry import Polygon
from shapely.geometry import Point

   
mpl.rc('font', family='sans-serif') 
mpl.rc('font',family='Times New Roman')
mpl.rc('text', usetex='false') 

fig = plt.figure(figsize=(8,8))


#########################################

ax = [plt.subplot(2,1,i+1) for i in range(2)]

for a in ax:
    #a.set_yticklabels(["netflux"])
   # a.set_ylabel('Netflux',fontsize=15)
    #a.set_xlabel('Energy', fontsize=15)
    #if a==ax[1]: #or a==ax[3]:
       # a.yaxis.tick_right()        
    for tick in a.xaxis.get_major_ticks():
        tick.label1.set_fontsize(18)
        tick.label1.set_fontweight('bold')#'bold
    for tick in a.yaxis.get_major_ticks():
        tick.label1.set_fontsize(18)             ############### Save thick axes
        tick.label2.set_fontsize(18)
        tick.label1.set_fontweight('bold')#
        tick.label2.set_fontweight('bold')#
    #a.set_yticklabels(np.arange(0,10,5), fontsize=15)
    
    #tick.label1.set_fontweight()
   # a.set_yticklabels(np.arange(0,10,5), fontsize=20)
    #a.set_ylim(0, 18)

plt.subplots_adjust(wspace=0, hspace=0)


#############
#%% 
############### check input EWs and Gammas #################


file=open("../ratios/fake_ratios/ew_fake.txt","r")
lines=file.readlines()
ew_org=[]
for x in lines:
    ew_org.append(float(x[:-1]))
file.close()


file=open("../ratios/fake_ratios/gamma_fake.txt","r")
lines=file.readlines()
gamma_org=[]
for x in lines:
    gamma_org.append(float(x[:-1]))
file.close()



############################################


#%%

xfull = fits.open('../catalog/CDWFS_I-3.6_v200226-final.fits')
xfull_dec = xfull[1].data['CDWFS_DEC']
xfull_ra = xfull[1].data['CDWFS_RA']
xfull_xid = xfull[1].data['CDWFS_ID']
xfull_z = xfull[1].data['zbest']
xfull_hbc = xfull[1].data['CDWFS_NET_H']

xfull_pos = np.vstack((xfull_xid, xfull_ra, xfull_dec, xfull_z, xfull_hbc)).T


#%% Calculate the fake array of crosses


ener = np.arange(0.075,15.0,0.15)
array_ratio=np.append(ener,ener)
array_ratio = np.append(array_ratio, [0]*xfull_pos.shape[1])
num=0


for j in range(0,len(xfull_pos)):

    xid = int(xfull_pos[j,0])
    z = xfull_pos[j,3]
    
    if os.path.isfile("../ratios/fake_ratios/"+str(xid)+"_fake.txt") == True: ## change directory ##
        
        
    
        f = open("../ratios/fake_ratios/"+str(xid)+"_fake.txt","r")
        lines=f.readlines()
        ener_obs = []
        flux_obs = []
        ener_obs_err = []
        flux_obs_err = []
        
        for x in lines[2:]:
            ener_obs.append(x.split(' ')[0])
            flux_obs.append(x.split(' ')[1])
            ener_obs_err.append(x.split(' ')[2])
            flux_obs_err.append(x.split(' ')[3][:-1])
            
            
        ener_obs = [float(i) for i in ener_obs] 
        flux_obs = [float(i) for i in flux_obs] 
        ener_obs_err = [float(i) for i in ener_obs_err] 
        flux_obs_err = [float(i) for i in flux_obs_err] 
              
        ener_res = np.asarray(ener_obs)*(1+z)
        flux_res = np.asarray(flux_obs)
        ener_res_err = np.asarray(ener_obs_err)*(1+z)
        flux_res_err = np.asarray(flux_obs_err)
        
        func_res = interp1d(ener_res, flux_res, kind='linear',bounds_error=False, fill_value=0)
        ener_rat = ener
        flux_rat = func_res(ener_rat)
       
        
        bsize_res = np.diff(ener_res)
        bsize_res = np.append(bsize_res,0)        
        fsize = interp1d(ener_res, bsize_res, kind='linear',bounds_error=False, fill_value=0)        
        bsize = fsize(ener)
        
        
        func_res_err = interp1d(ener_res, flux_res_err, kind='linear',
                                bounds_error=False, fill_value=0)
        
        err = func_res_err(ener)
        
        flux_rat_err = err*np.sqrt(bsize/0.15)

        
            
        #ener_tot = cbind(ener_rat,ener_rat_err)
        flux_rat = np.append(flux_rat, flux_rat_err)
        

        flux_rat = np.append(flux_rat, xfull_pos[j,:])
        

        array_ratio = np.vstack((array_ratio, flux_rat))
        num=num+1
        
    #else:
 
    #    print(j)
        
    j=j+1



cdwfs_fake = array_ratio[1:,:] 

df_fake = pd.DataFrame.from_records(cdwfs_fake)

df_fake.rename(columns={200: 'xid', 201:'xra',202:'xdec',203:'z', 204:'hbc'}, inplace= True)

df_fake['fake_ew'] =  ew_org
df_fake['fake_gamma'] =  gamma_org

#%% Cross Matched Sources

def xidcrossmatch(fir_num):
    crossmatch = np.load("../array_"+str(fir_num)+"_ratio.npy")
    crossmatch = crossmatch[1:,200:206]
    return crossmatch



cdwfs_70 = xidcrossmatch(70)
cdwfs_100 = xidcrossmatch(100)
cdwfs_160 = xidcrossmatch(160)
cdwfs_250 = xidcrossmatch(250)
cdwfs_350 = xidcrossmatch(350)
cdwfs_500 = xidcrossmatch(500)

cdwfs_cross1 = np.vstack((cdwfs_70,cdwfs_100))
cdwfs_cross2 = np.vstack((cdwfs_160,cdwfs_250))
cdwfs_cross3 = np.vstack((cdwfs_350,cdwfs_500))
cdwfs_cross4 = np.vstack((cdwfs_cross1,cdwfs_cross2))

cdwfs_cross = np.vstack((cdwfs_cross4,cdwfs_cross3))

cdwfs_uni = cdwfs_cross[0:2,:]
for i in range(0,len(cdwfs_cross)):
    xid = cdwfs_cross[i,5]
    if (xid in cdwfs_uni[:,5]) == False:
        cdwfs_uni = np.vstack((cdwfs_uni,cdwfs_cross[i,:]))
    i=i+1

df_cross = pd.DataFrame.from_records(cdwfs_uni)
df_cross.rename(columns={0:'hbc', 1:'z',  2:'fir', 3:'xra', 4:'xdec', 5:'xid'},
          inplace=True)


# Calculate LIR
ldis  = cosmo.luminosity_distance(df_cross["z"]).value # in Mpc
df_cross["lir"]  = 4*np.pi*ldis**2*df_cross["fir"] # Unit is not standard

df_cross = df_cross.drop('hbc', axis=1)

#df_fake = pd.concat([df_ratio, df_xid], axis=1) #side by side


#%% Creating two weights

ener_ind = np.arange(0,100,1)
  
fake_ener = df_fake[ener_ind]

# for i in range(0,df_fake.shape[0]):

  
#   inte = integrate.quad(f,0,val)
#   f = interp1d(ener, ratioi )
  
#   inte = integrate(approxfun(ener, ratioi,method="linear"), lower = ener[26], upper = ener[53])$value
  
  
#   inte_all <- c(inte_all, inte)
  
# }

# inte_log <- log10(abs(inte_all))
# norm_func <- inte_all/(inte_log-0.99*min(inte_log[is.finite(inte_log)]))^3
# weight1 = (inte_log-0.99*min(inte_log[is.finite(inte_log)]))^3

weight2 = ener/ener


#%% Judge by field

all_pts = np.asarray(df_fake[["xra","xdec"]])
cross_pts = np.asarray(df_cross[["xra","xdec"]])


#sp_ot = pd.read_csv('../figures/spire_outer.csv')
sp_ot = pd.read_csv('../figures/spire_inner.csv')
poly_ot = Polygon(np.asarray(sp_ot))



def inorout(points_arr, poly):
    in_ind=[]
    #out_ind=[]
    for i in range(0,len(points_arr)):
        x = points_arr[i,]
        if Point(x).within(poly) == True:
            in_ind.append(i)
            
        #else:
        #    out_ind.append(i)
        
        i=i+1
    return in_ind#, out_ind


sp_ot_cross_ind = inorout(cross_pts, poly_ot)
sp_ot_all_ind = inorout(all_pts, poly_ot)


#%% Detected vs Undetected -------------------------------------------


ot_detect = df_cross.loc[sp_ot_cross_ind]
ot_all = df_fake.loc[sp_ot_all_ind]


ot_all= pd.merge(ot_all, ot_detect, on='xid', how='outer')
ot_all.rename(columns={'xra_x': 'xra', 'xdec_x':'xdec','z_x':'z'}, inplace= True)
ot_all.drop(columns = {'xra_y','xdec_y','z_y'}, inplace=True)


ot_detect = ot_all.dropna(subset=['lir'])
ot_undetect = ot_all[~ot_all.xid.isin(ot_detect.xid)]

############## Gamma Clip #################

# gammacut = 1.8
# xid_replace = ot_all[ot_all['fake_gamma'] > gammacut]['xid']

# with open("../ratios/fake_ratios/xid_replace.txt", "w") as f:
#     for s in xid_replace:
#         s = int(s)
#         f.write(str(s) +"\n")      
# f.close()



#%% Cuts and Selections -------------------------------------------

# redshift cut
ot_detect = ot_detect[0.5 < ot_detect['z']]
ot_detect = ot_detect[ot_detect['z'] < 1.5]
ot_undetect = ot_undetect[0.5 < ot_undetect['z']]
ot_undetect = ot_undetect[ot_undetect['z']<1.5]
print('There is a redshift cut')
#print('There is no redshift cut')

# count cut
cts_cut = 30
ot_detect = ot_detect[ot_detect['hbc']>cts_cut]
ot_undetect = ot_undetect[ot_undetect['hbc']>cts_cut]

print('The hbc cut is %i' % cts_cut)


print('%i undetected sources and %i detected sources'
      % (ot_undetect.shape[0],ot_detect.shape[0] ) )

print('Input EW mean of undetected source is %.2f'
      % (ot_undetect['fake_ew'].mean() ))

print('Input EW mean of detected source is %.2f'
      % (ot_detect['fake_ew'].mean() ))


#%% Cuts and Selections -------------------------------------------

dt_ave = ot_detect.loc[:,range(26,55)].mean(axis = 1)
un_ave = ot_undetect.loc[:,range(26,55)].mean(axis = 1)

dt_filter = 1-sigma_clip(dt_ave,sigma=3).mask
un_filter = 1-sigma_clip(un_ave, sigma=3).mask

ot_detect['clip'] = dt_filter
ot_undetect['clip'] = un_filter

ot_detect = ot_detect[ot_detect['clip'] ==1]
ot_undetect = ot_undetect[ot_undetect['clip'] ==1]

# fig=plt.figure(figsize=(8,5))

# for i in range(0, ot_detect.shape[0]):
#     plt.plot(ener, ot_detect.iloc[i,0:100],'.')

# plt.show()



#%% Stacking -------------------------------------------


flux = (ot_undetect[ener_ind] * weight2)
f_ot_low = flux.sum(axis = 0, skipna = True)/flux.shape[0]




flux = (ot_detect[ener_ind] * weight2)
f_ot_high = flux.sum(axis = 0, skipna = True)/flux.shape[0]

                
#%%  Fit  

def fefit(ydata):

    index=np.array([6,7,8,9,10,11,12,13,14,15,16])

    xdata = ener[32:55]
    ydata = ydata[32:55]
    y_noise = ydata*0.1
    
    new_x = np.delete(xdata, index)
    new_y = np.delete(ydata, index)
    new_ynoise = np.delete(y_noise, index)
    
    def func2(x, a, b):
        return a*x+b


    p=np.array([0.5,0.3])
    
    popt_pl, pcov = curve_fit(func2, new_x, new_y, sigma=new_ynoise, p0=p)


    def func(x,c, d):
        para = popt_pl[0]
        parb = popt_pl[1]

        return c*np.exp(-0.5*((x-6.4)/d)**2)+para*x+parb


    popt, pcov = curve_fit(func, xdata, ydata, sigma=y_noise, p0=p,bounds=((0, 0), (np.inf, 0.5)))


    fungau = popt[0]*np.exp(-0.5*((xdata-6.4)/(popt[1]))**2)
    funpl = popt_pl[0]*xdata+popt_pl[1]
    funmod = fungau+funpl

#def integrand(x):
#    return popt[0]*np.exp(-0.5*((x-6.4)/popt[1])**2)
    def integrand(x):
        return popt[0]*np.exp(-0.5*((x-6.4)/popt[1])**2) / (popt_pl[0]*x+popt_pl[1])
    
    ans, err = quad(integrand, 5., 8.)
    y_line=popt_pl[0]*6.4+popt_pl[1]
    
    width = ans/y_line   

    
    return xdata, ydata, fungau, funmod, popt, width


xdata, ydata, fungau, funmod, popt, width = fefit(np.asarray(f_ot_low))


ax[0].plot(xdata, fungau, "b--", linewidth=1)
ax[0].plot(xdata, funmod, "-", color="purple",linewidth=4,label="Best Fits")
ax[0].errorbar(xdata, ydata, yerr=ydata*0.1,fmt='.', ecolor="blue",capsize=4, markersize=10,label="Net Flux")

ax[0].set_xlim(5.01, 7.99)

print('EW and Sigma of undetected sources is %.2f and %.2f' 
      % (width, popt[1]))

# #############################################################################################################################################################################################################

xdata, ydata, fungau, funmod, popt, width = fefit(np.asarray(f_ot_high))

ax[1].plot(xdata, fungau, "b--", linewidth=1)
ax[1].plot(xdata, funmod, "-", color="purple",linewidth=4,label="Best Fits")
ax[1].errorbar(xdata, ydata, yerr=ydata*0.1,fmt='.', ecolor="blue",capsize=4, markersize=10,label="Net Flux")


ax[1].set_xlim(5.01, 7.99)
print('EW and Sigma of detected sources is %.2f and %.2f' 
      % (width, popt[1]))


plt.show()


#%% try different seperation

try_low = ot_detect[ot_detect['lir'] < ot_detect['lir'].mean()]
try_high = ot_detect[ot_detect['lir'] >= ot_detect['lir'].mean()]

flux = (try_low[ener_ind] * weight2)
f_ot_low = flux.sum(axis = 0, skipna = True)/flux.shape[0]


flux = (try_high[ener_ind] * weight2)
f_ot_high = flux.sum(axis = 0, skipna = True)/flux.shape[0]





#%%  Fit 
# ##################
# print(np.mean(ewl))
# print(np.mean(ewh))
# #######################



# pl.show()


# pp = PdfPages('fake.pdf')
# pp.savefig(fig)
# pp.close()