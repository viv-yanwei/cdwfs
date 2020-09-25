##!/bin/sh
#  
#
##!/bin/sh
#  
#import numpy as np
#
#import os
#import os.path
#from astropy.io import fits
#import math
#
#
#
#
#file=open("selected_xid.txt","r")
#lines=file.readlines()
#xid=[]
#for x in lines:
#    xid.append(x[:-1])
#file.close()
#
#
#
#file=open("selected_xid.txt","r")
#lines=file.readlines()
#xid=[]
#for x in lines:
#    xid.append(x[:-1])
#file.close()

# 2–7 keV (hard band; HB)

#gal = fits.open("/Users/yan/Box/CDFS/HostGalaxy/7ms_her_s15_1903.fits")
#
#hbf = gal[1].data["HARD_BAND_FLUX"]
#
#xid_7ms = gal[1].data["XID_SOURCE_NUMBER"]
#xid_7ms=list(xid_7ms)
#
#
#for i in range(0,len(xid)):
#    
#    os.chdir('/Users/yan/Box/CDFS/1119/source/'+xid[i])
#
#    ind = [xid_7ms.index(t) for t in xid_7ms if t == float(xid[0])]
#    
#    with open('hbf.txt', 'a') as f:
#        f.write(str(hbf[ind[0]]) +"\n")
#    
#    i=i+1
#
#


import numpy as np
from astropy.io import fits
import os
import os.path
import math
import random as rd

from ciao_contrib.runtool import *
from sherpa.astro.ui import *  

full = fits.open("/Users/yan/Dropbox/sharinig_research/CDWFS/catalog/CDWFS_I-3.6_v200226-final.fits")
#
cid = full[1].data["CDWFS_ID"]
zbest = full[1].data["zbest"]
fflux =  full[1].data["CDWFS_FLUX_F"]

czf = zip(cid,zbest,fflux)

czf_z=[] # 1st cut - redshift: 
         # select those with normal redshift (aka z>0)
for x in czf:
    if x[1]>0:
        czf_z.append(x)   # 3 col: xid, z, full_band_flux



nodata=[]
ew_fake=[]
gamma_list=[]

for i in range(0, len(czf_z)):
    xid = str(czf_z[i][0])
    
    if os.path.isfile("../multiple_obs/"+xid+"_combined_src.pi") == True:

        arfname = "../multiple_obs/"+xid+"_combined_src.arf"
        rmfname = "../multiple_obs/"+xid+"_combined_src.rmf"
        srcname = "../multiple_obs/"+xid+"_combined_src.pi"

        z=czf_z[i][1]
        fflux = czf_z[i][2]

    
    elif os.path.isfile("../single_obs/"+xid+".pi") == True:

        arfname = "../single_obs/"+xid+".arf"
        rmfname = "../single_obs/"+xid+".rmf"
        srcname = "../single_obs/"+xid+".pi"

        z=czf_z[i][1]
        fflux = czf_z[i][2]

    


    else:
        nodata.append(xid)
        #print(xid)
        continue
    

    obs_flux = fflux
 

    set_source(1, const1d.ct*(powlaw1d.p*xszphabs.q+gauss1d.line))
    
    ct.c0 =  1
    
   
    p.gamma = 1.4
    
    
    q.nH=1
    q.redshift=z
    
   
    line.pos.min = 0
    line.pos.max = 100
    
    line.pos=6.4/(1+z)

    
    line.pos.min = 6.3/(1+z)
    line.pos.max = 6.5/(1+z)
     
    
    line.fwhm = 0.8/(1+z)  # change value 0.8 based on broadening
    line.ampl = 1
    


#    bkg1 = unpack_pha("1119_combined_bkg.pi")
    
        
    load_pha(2,srcname)
    t_expo = get_exposure(2)
    

    

    fake_pha(1,arfname, rmfname, t_expo)
#    fake_pha(1,arfname, rmfname, t_expo, bkg=bkg1)
#    subtract(1)


    ew0=eqwidth(p*q, p*q+line,id=1)

    gmean = 0.5
    ew_set = rd.gauss(gmean,0.2)
    if ew_set < 0 :
        ew_set=0
    
    
    ew_fake.append(ew_set)
    line.ampl = (ew_set/(1+z))/ew0
    
    
    #print(eqwidth(p*q, p*q+line,id=1))
    
    
    mdl_flux = calc_energy_flux(2,7,id=1)
    
    ct.c0 = obs_flux/mdl_flux

#    fake_pha(1,arfname, rmfname, t_expo, bkg=bkg1)
    fake_pha(1,arfname, rmfname, t_expo)

    
    

########## BKG #############

	#load_bkg("1119_combined_bkg.pi")
	#if os.path.isfile("spec_combined_bkg.arf") == True:
    #    	load_bkg_arf("spec_combined_bkg.arf")
	#else:
    #    	redo.append(xid)
	#if os.path.isfile("spec_combined_bkg.rmf") == True:
    #    	load_bkg_rmf("spec_combined_bkg.rmf")
	#else:
    #    	redo.append(xid)

#    subtract(1)

    group_counts(1,5)
    
#    create_model_component("powlaw1d","mdl1")
#    create_model_component("xsphabs", "mdl2")

    
    notice(1,10)
    ignore(6/(1+z),7/(1+z))
    
 
    set_source(1,powlaw1d.p*xszphabs.q+gauss1d.line)
 
    thaw(line.pos)
    thaw(line.fwhm)
    thaw(line.ampl)
    thaw(q.nH)
    thaw(p.gamma)
    thaw(p.ampl)
    
 
    line.pos=6.4/(1+z)
    line.fwhm = 0.1
    line.ampl=0
    
    
    freeze(line.ampl)
    freeze(line.pos)
    freeze(line.fwhm)

  
    p.gamma = 1.4
    
    q.redshift=z


    fit()
    plot_fit()
    notice(1,10)
    plot_data()
    r=get_ratio_plot()
    gamma_list.append(p.gamma.val)
    #print('Gamma parameter is: ', p.gamma.val)
    
    
    freeze(p.gamma)
    freeze(p.ampl)
    freeze(q.nH)
    
    thaw(line.ampl)
    
    fit()
    
    #print(eqwidth(p*q, p*q+line,id=1))
    #print(i)
    
    
    save_arrays("../ratios/fake_ratios/"+xid+"_fake.txt", [r.x,r.y,r.xerr,r.yerr], ["x","y","xerr","yerr"],clobber="TRUE")
        
        
    i=i+1

    


########## files with problems


with open("../ratios/fake_ratios/nodata.txt", "w") as f:
    for s in nodata:
    		f.write(str(s) +"\n")
            
f.close()

with open("../ratios/fake_ratios/ew_fake.txt", "w") as f:
    for s in ew_fake:
        f.write(str(s) +"\n")

            
f.close()

with open("../ratios/fake_ratios/gamma_fake.txt", "w") as f:
    for s in gamma_list:
    		f.write(str(s) +"\n")
            
f.close()
