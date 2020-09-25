setwd('/Users/yan/Dropbox/sharinig_research/CDWFS/codes')
library(reticulate)
library(rlist)
library(FITSio)
library(magicaxis)
library(ggplot2)
library(astrolibR)
library(gtable)
library(grid)
library(dplyr)
library(sp)

np <- import("numpy")
#cdwfs <- np$load("../array_ratio.npy") # full
cdwfs_70 <- np$load("../array_70_ratio.npy")
ener <- cdwfs_70[1,1:100]
cdwfs_70 <- cdwfs_70 [-1,]
cdwfs_100 <- np$load("../array_100_ratio.npy")
cdwfs_100  <- cdwfs_100 [-1,]
cdwfs_160 <- np$load("../array_160_ratio.npy")
cdwfs_160 <- cdwfs_160 [-1,]
cdwfs_250 <- np$load("../array_250_ratio.npy")
cdwfs_250 <- cdwfs_250 [-1,]
cdwfs_350 <- np$load("../array_350_ratio.npy")
cdwfs_350 <- cdwfs_350 [-1,]
cdwfs_500 <- np$load("../array_500_ratio.npy")
cdwfs_500 <- cdwfs_500 [-1,]


cdwfs <- rbind(cdwfs_70,cdwfs_100,cdwfs_160,cdwfs_250,cdwfs_350,cdwfs_500)
# [1:100] ratio
# [101:200] err
# [201:206] hbc,z,fir,xra,xdec,xid

# Remove duplicate rows of the dataframe
uniq <- cdwfs[1:2,]
for (i in c(1:nrow(cdwfs))) {
  if (cdwfs[i,206] %in% uniq[,206] == FALSE) {
    uniq <- rbind(uniq,cdwfs[i,])
  }
  i=i+1
}

cdwfs <- uniq

sp_ot <- read.csv(file = '../figures/spire_outer.csv')
sp_ot <- data.matrix(sp_ot, rownames.force = NA)

sp_in <- read.csv(file = '../figures/spire_inner.csv')
sp_in <- data.matrix(sp_in, rownames.force = NA)



# Creating two weights ----------------------------------------------------

inte_all <- NULL
for (i in c(1:nrow(cdwfs))) {
  
  ratioi <- cdwfs[i,1:100]
  inte <- integrate(approxfun(ener, ratioi,method="linear"), lower = ener[26], upper = ener[53])$value
  
  
  inte_all <- c(inte_all, inte)
  
}

inte_log <- log10(abs(inte_all))
norm_func <- inte_all/(inte_log-0.99*min(inte_log[is.finite(inte_log)]))^3
weight1 <- (inte_log-0.99*min(inte_log[is.finite(inte_log)]))^3

weight2 <- rep(1,nrow(cdwfs))


# Clips -------------------------------------------------------------------

# Calculate LIR first
ldis <-lumdist(cdwfs[,202]) # in Mpc
lir <- 4*pi*ldis^2*cdwfs[,203] # Unit is not standard
#
#example: spec_z <- spec[ which(spec[,107] < 1.5 & spec[,107] > 0.5 ), ]

cdwfs <- cbind(cdwfs, lir)  # 207 lir


# Calculate location
#0: point is strictly exterior to pol; 
#1: point is strictly interior to pol; 
#2: point lies on the relative interior of an edge of pol; 
#3: point is a vertex of pol.
loc_outer <- point.in.polygon(nf_rat[,204], nf_rat[,205], 
                              sp_ot[,1], sp_ot[,2], mode.checked=FALSE)

loc_inner <-point.in.polygon(cdwfs[,204], cdwfs[,205], 
                             sp_in[,1], sp_in[,2], mode.checked=FALSE)

cdwfs <- cbind(cdwfs, loc_inner)
cdwfs <- cbind(cdwfs, loc_outer)  # add 208,209 

#ave_hb<- apply(cdwfs_hc[,45:65], 1, mean)
#cdwfs_clip <- cdwfs_hc[which(magclip(ave_hb, sigma = 3)$clip==TRUE),]



# Choose Weight and Calculation -------------------------------------------


cdwfs_loc <- cdwfs[which(cdwfs[,208]==1),] # inner region
 
cdwfs_loc <- cdwfs[which(cdwfs[,209]==1),] # outer region

weight <- weight2
nf_rat <- NULL

for  (i in c(1:nrow(cdwfs_loc))) {
  nfi <- cdwfs_loc[i,1:100] * weight[i]
  nferri= ( cdwfs_loc[i,101:200] * weight[i])^2
  
  nf_rati <- append(nfi, nferri)
  nf_rati <- append(nf_rati, cdwfs_loc[i,201:209])   # add props

  nf_rat <- rbind(nf_rat, nf_rati)
}

# Subgroups ---------------------------------------------------------------

# select from nf_rat
# 1. outer detected
cdwfs_z <- nf_rat[which(nf_rat[,202]<=2),]
cdwfs_hc <- cdwfs_z[which(cdwfs_z[,201]>80),]

cdwfs_out_d <- cdwfs_hc


# 2. outer non-detected
uniq <- rep.int(1,209)
for (i in c(1:nrow(nf_rat))) {
  if (nf_rat[i,206] %in% cdwfs_hc[,206] == FALSE) {
    
    uniq <- rbind(uniq,nf_rat[i,])
  }
  i=i+1
}

cdwfs_out_nd <- uniq[-1,] 



# full
# nf_clip <- colSums (nf_rat[,1:100], na.rm = TRUE, dims = 1) / nrow(nf_rat)
# clip <- data.frame(x=ener, y=nf_clip)
# 
# 
# nf_clip_err <- sqrt(colSums (nf_rat[,101:200], na.rm = TRUE, dims = 1))/nrow(nf_rat)
# clip_err <- data.frame(x=ener, upper=nf_clip+nf_clip_err/2, lower=nf_clip-nf_clip_err/2)


######### Plot ####


#Plot full


# g_clip <- ggplot(clip, aes(x=ener, y=nf_clip)) 
# g_clip <- g_clip + xlab("Energy (keV)") + ylab("Flux")
# g_clip <- g_clip + geom_point(size=1, colour="blue")
# g_clip <- g_clip + geom_errorbar(clip_err, mapping=aes(x=ener, ymin=upper, ymax=lower), width=0.2, size=0.5, color="blue")
# g_clip <- g_clip+ xlim(5,8) #+ ylim(0.,)
# g_clip
# 
# 
# f_sel <- cbind(ener, nf_clip, nf_clip_err)
# writeFITSim(f_sel, file = "f_sel.fits", type = "double", bscale = 1, bzero = 0, c1 = NA, c2 = NA,
#            crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
#            axDat = NA, header = '')
# 


# seperate by nondetections:

irl <- colSums (cdwfs_out_nd[,1:100], na.rm = TRUE, dims = 1) / nrow(cdwfs_hc)
irh <- colSums (cdwfs_out_d[,1:100], na.rm = TRUE, dims = 1) / nrow(cdwfs_out_nd)

df_irl <- data.frame(x=ener, y=irl)
df_irh <- data.frame(x=ener, y=irh)


# seperate by IR (,204)

#nf_rat <- nf_rat[which(nf_rat[,202] < 1.5 & nf_rat[,202] > 0.5 ),]
# mlir <- median(nf_rat[,204])
# nf_irl <- nf_rat[which(nf_rat[,204]<mlir),]
# nf_irh <- nf_rat[which(nf_rat[,204]>=mlir),]
# 
# irl <- colSums (nf_irl[,1:100], na.rm = TRUE, dims = 1) / nrow(nf_irl)
# irh <- colSums (nf_irh[,1:100], na.rm = TRUE, dims = 1) / nrow(nf_irh)
# 
# df_irl <- data.frame(x=ener, y=irl)
# df_irh <- data.frame(x=ener, y=irh)
#nf_clip_err <- sqrt(colSums (nf_rat[,101:200], na.rm = TRUE, dims = 1))/nrow(nf_rat)
#clip_err <- data.frame(x=ener, upper=nf_clip+nf_clip_err/2, lower=nf_clip-nf_clip_err/2)


#Plot ir


g_irl <- ggplot(df_irl, aes(x=ener, y=irl))
g_irl <- g_irl + xlab("Energy (keV)") + ylab("Flux")
g_irl <- g_irl + geom_point(size=1, colour="blue")
#g_irl <- g_clip + geom_errorbar(clip_err, mapping=aes(x=ener, ymin=upper, ymax=lower), width=0.2, size=0.5, color="blue")
g_irl <- g_irl+ xlim(2,10) #+ ylim(0.,)
g_irl

g_irh <- ggplot(df_irh, aes(x=ener, y=irh))
g_irh <- g_irh + xlab("Energy (keV)") + ylab("Flux")
g_irh <- g_irh + geom_point(size=1, colour="blue")
#g_irl <- g_clip + geom_errorbar(clip_err, mapping=aes(x=ener, ymin=upper, ymax=lower), width=0.2, size=0.5, color="blue")
g_irh <- g_irh+ xlim(2,10) #+ ylim(0.,)
g_irh


f_irl <- cbind(ener, irl)
writeFITSim(f_irl, file = "f_irl.fits", type = "double", bscale = 1, bzero = 0, c1 = NA, c2 = NA,
           crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
           axDat = NA, header = '')

f_irh <- cbind(ener, irh)
writeFITSim(f_irh, file = "f_irh.fits", type = "double", bscale = 1, bzero = 0, c1 = NA, c2 = NA,
            crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
            axDat = NA, header = '')








######### Bootstrapping ####
# 
# boot=NULL
# for (i in c(1:500)) {
#   index <- sample(1:nrow(nf_rat[,1:100]),replace = TRUE)
#   resample=NULL
#   for (j in c(1:length(index))){
#     resample <- rbind(resample, nf_rat[index[j],1:100])
#   }
#   nf_resample <- colSums(resample, na.rm = TRUE, dims = 1) / nrow(resamples)
#   boot <- cbind(boot, nf_resample)
# }
# 
# 
# writeFITSim(boot, file = "f_boot.fits", type = "double", bscale = 1, bzero = 0, c1 = NA, c2 = NA,
#             crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
#             axDat = NA, header = '')


boot_irl=NULL
for (i in c(1:500)) {
  index <- sample(1:nrow(nf_irl[,1:100]),replace = TRUE)
  resample=NULL
  for (j in c(1:length(index))){
    resample <- rbind(resample, nf_irl[index[j],1:100])
  }
  nf_irl_resample <- colSums(resample, na.rm = TRUE, dims = 1) / nrow(resample)
  boot_irl <- cbind(boot_irl, nf_irl_resample)
}


writeFITSim(boot_irl, file = "f_irl_boot.fits", type = "double", bscale = 1, bzero = 0, c1 = NA, c2 = NA,
            crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
            axDat = NA, header = '')


boot_irh=NULL
for (i in c(1:500)) {
  index <- sample(1:nrow(nf_irh[,1:100]),replace = TRUE)
  resample=NULL
  for (j in c(1:length(index))){
    resample <- rbind(resample, nf_irh[index[j],1:100])
  }
  nf_irh_resample <- colSums(resample, na.rm = TRUE, dims = 1) / nrow(resample)
  boot_irh <- cbind(boot_irh, nf_irh_resample)
}


writeFITSim(boot_irh, file = "f_irh_boot.fits", type = "double", bscale = 1, bzero = 0, c1 = NA, c2 = NA,
            crpixn = NA, crvaln = NA, cdeltn = NA, ctypen = NA, cunitn = NA,
            axDat = NA, header = '')

