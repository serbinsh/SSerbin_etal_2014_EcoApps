####################################################################################################
# Supplement 2
#
# Shawn P. Serbin, Aditya Singh, Brenden E. McNeil, Clayton C. Kingdon, and Philip A. Townsend. 2014.
# Spectroscopic determination of leaf morphological and biochemical traits for northern temperate and 
# boreal tree species. Ecological Applications.
#            
#
# Simplified script file illustrating the methods for calibrating the dry spectr PLSR models and 
# developing the uncertainty estimates.
#
#
# Programmer: Shawn P. Serbin, October, 2013
# R version: Up to 3.03
####################################################################################################


#--------------------------------------------------------------------------------------------------#
# Close all devices and delete all variables.
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- Load required libraries ---------------------------------------------------------#
# Info: Loads required R libraries and warns if package is not availible.
ok = require(pls) ; if (! ok) 
  stop("*** Package pls is not available.  This is needed for model optimization ***")

ok = require(plotrix) ; if (! ok) 
  stop("*** Package plotrix is not available.  This is needed visualization of results ***")

# Script options
pls.options(plsralg = "oscorespls")
pls.options("plsralg")
rm(ok)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Set ouput directory
out.dir = '/Users/serbin/DATA/Dropbox/MANUSCRIPTS/FERST_LAB/Dissertation_Manuscripts/FFT_Dry_Spectra_Chem_Paper/R_output/Leaf_LMA/'
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Import dry spectra dataset
spec.dir <- '/Users/serbin/DATA/Dropbox/MANUSCRIPTS/FERST_LAB/Dissertation_Manuscripts/FFT_Dry_Spectra_Chem_Paper/Data/'
dry.spectra <- read.table(paste(spec.dir,'NASA_FFT_DS_Refl_Spectra_v2.csv',sep=""), header=TRUE,sep=",")
dry.spectra[1:30,1:10]
rm(spec.dir)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Generate diagnostic spectra figure?

Create.fig <- FALSE

if (Create.fig) {
  ### Create dry spec diagnostic figure
  dims <- dim(dry.spectra)
  spec <- dry.spectra[,3:dims[2]]
  mean.spec <- colMeans(spec)
  spec.quant <- apply(spec,2,quantile,probs=c(0.05,0.95))
  
  # View/Plot spectra
  pdf(paste(out.dir,'FFT_Raw_Dry_Spectra_QC.pdf',sep=""),height=7,width=8)
  par(mar=c(4,4,1,1.2)) #B,L,T,R
  matplot(t(spec), type = "l", lty = 1, ylab = "Reflectance (%)", xaxt = "n")
  ind <- pretty(seq(from = 350, to = 2500, by = 1)) # Using pretty to standardize the axis
  ind <- ind[ind >= 350 & ind <= 2500]
  ind <- (ind - 349) / 1
  axis(1, ind, colnames(spec)[ind]) # add column names to wavelengths
  # Mean spectra
  lines(mean.spec,lwd=9)
  # lines(mean.spec+(sd.spec*1.96),lty=3,lwd=5,col="black")
  # lines(mean.spec-(sd.spec*1.96),lty=3,lwd=5,col="black")
  # CIs
  lines(spec.quant[1,],lty=1,lwd=7,col="dark grey")
  lines(spec.quant[2,],lty=1,lwd=7,col="dark grey")
  legend("bottomright",legend=c("Mean","95% CI"),lty=c(1,1),
         col=c("black","dark grey"),lwd=3)
  box(lwd=2.2)
  dev.off()
  rm(spec,mean.spec,spec.quant,ind)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Import LMA data
lma.dir = '/Users/serbin/DATA/Dropbox/MANUSCRIPTS/FERST_LAB/Dissertation_Manuscripts/Data/SLA_LMA/'
lma.data = read.table(paste(lma.dir,'FFT_SLA_LMA_DATA_v2_4R_Updated.csv',sep=""),header=TRUE,sep=",")
names(lma.data)
dims <- dim(lma.data)
dims

### Sites
unique(lma.data$SITE)

### Species
unique(lma.data$SPECIES)

hist(lma.data$LMA_g_DW_m2,breaks=20)
#--------------------------------------------------------------------------------------------------#


#---------------- Merge data with spectra ---------------------------------------------------------#
merged.cal <- merge(lma.data, dry.spectra,by.x=c("SAMPLE_NAME","SAMPLE_YEAR"),
                   by.y=c("Sample_Name","Sample_Year"),
                   all.x=FALSE,all.y=FALSE)
dim(merged.cal)
# 766, 2167

hist(merged.cal$LMA_g_DW_m2)

unique(merged.cal$SITE)
unique(merged.cal$SPECIES)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Spectral/chemical correlations
waves <- data.frame(seq(1,2151,1),seq(350,2500,1))
dims <- dim(merged.cal)

# Raw spectra
spec <- as.matrix(merged.cal[,17:dims[2]])
# Raw
spec_corr = cor(as.matrix(merged.cal[,17:dims[2]]), merged.cal$LMA_g_DW_m2)
#Log 1/R transform
#spec_corr = cor(as.matrix(log10(1/merged.cal[,17:dims[2]])), merged.cal$LMA_g_DW_m2)

pdf(paste(out.dir,'FFT_Leaf_LMA_Spectra_Correlations.pdf',sep=""),height=12,width=8)
par(mfrow=c(2,1),mar=c(4,4.6,1,1.4)) #B, L, T, R
matplot(t(spec), type = "l", lty = 1, ylab = "Reflectance (%)", xaxt = "n",ylim=c(0,0.9))
ind <- pretty(seq(from = 350, to = 2500, by = 1)) # Using pretty to standardize the axis
ind <- ind[ind >= 350 & ind <= 2500]
ind <- (ind - 349) / 1
axis(1, ind, colnames(spec)[ind]) # add column names to wavelengths

plot(waves[,2],spec_corr,xlab="WAVELENGTH (nm)",ylab="CORRELATION",
     main="Leaf LMA (gDW/m2)", cex=0.01)
lines(waves[,2],spec_corr,lwd=4)
abline(h=0,lty=2,lwd=1.5,col="grey80")
box(lwd=2)
dev.off()

# Output correlation data
spec_corr <- data.frame(spec_corr)
names(spec_corr) <- c("Correlation")
write.csv(spec_corr,paste(out.dir,'FFT_Leaf_LMA_Spectra_Correlations.csv',sep=""),
          row.names=TRUE)
rm(spec,spec_corr,waves,ind)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Subset spectra for analysis 
dims <- dim(merged.cal)
temp <- merged.cal[,17:dims[2]]
dry.spectra2 <- as.matrix(temp[,151:2051]) #151:2101 (500-2400nm), 551:2101 (900 - 2400nm),
# 851:2051 (1200 - 2400nm)
dry.spectra2[1:5,1:10]
dry.spectra2[1:5,1899:1901]

rm(temp, dry.spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Build full PLSR dataset 
full.plsr.data <- data.frame(merged.cal[,1:16],dry.spectra2)
full.plsr.data[1:5,1:16]
dim(full.plsr.data)
rm(merged.cal,dry.spectra2)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Subset data into cal/val by site
sites <- unique(full.plsr.data$SITE)
create.seed <- TRUE  #TRUE/FALSE

# Sample proportion for cal data
prop <- 0.8

# Random seed
if (create.seed){
  set.seed(as.vector(round(runif(5,min=5,max=9))))
  ### Write out seed
  .Random.seed[1:6]
  seed.save <- .Random.seed
  write.table(seed.save,paste(out.dir,"random.seed",sep=""));
}

### Read in previous random seed
seed <- read.table(paste(out.dir,"random.seed",sep=""))[,1];
.Random.seed <- seed

cal.plsr.data <- 0
val.plsr.data <- 0
j <- 1
for (i in sites){
  print(paste("Site: ",i,sep=""))
  temp.data <- full.plsr.data[which(full.plsr.data$SITE==i),]
  rows <- sample(1:nrow(temp.data),floor(prop*nrow(temp.data)))
  cal_data = droplevels(temp.data[rows,])
  val_data = droplevels(temp.data[-rows,])
  
  if(j==1){
    cal.plsr.data <- cal_data
    val.plsr.data <- val_data
  } else {
    cal.plsr.data <- rbind(cal.plsr.data,cal_data)
    val.plsr.data <- rbind(val.plsr.data,val_data)
  }
  
  j <- j+1
}
rm(temp.data)

# Datasets:
# cal.plsr.data -- For building PLSR model
print(paste("Cal observations: ",dim(cal.plsr.data)[1],sep=""))
# val.plsr.data -- Independent (external) PLSR model validation data ~20 of data
print(paste("Val observations: ",dim(val.plsr.data)[1],sep=""))

pdf(paste(out.dir,'FFT_Leaf_LMA_Cal_Val_Histograms.pdf',sep=""),height=12,width=8)
par(mfrow=c(2,1),mar=c(4,4.6,1,1.4)) #B, L, T, R
hist(cal.plsr.data$LMA_g_DW_m2)
hist(val.plsr.data$LMA_g_DW_m2)
dev.off()

write.csv(full.plsr.data,file=paste(out.dir,"FFT_Leaf_LMA_Full_PLSR_Dataset.csv",sep=""),row.names=FALSE)

rm(cal_data,val_data,i,j,prop,rows,sites)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Run calibration PLSR analysis to select optimal number of components
cal.plsr.data <- cal.plsr.data2
rm(cal.plsr.data2,outlier)
dims = dim(cal.plsr.data)

k <- round(dims[1]/10)
segs = cvsegments(dims[1],k = k, type="random")
LeafLMA.pls = plsr(LMA_g_DW_m2~Spectra,scale=FALSE,ncomp=15,validation="CV",segments=segs,
                   trace=TRUE, method = "oscorespls", data=cal.plsr.data)
rm(segs)

# Examine raw PLSR output
summary(LeafLMA.pls)
plot(RMSEP(LeafLMA.pls), legendpos = "topright")
names(LeafLMA.pls)

predplot(LeafLMA.pls, ncomp = 8:14, asp = 1, line = TRUE,which = c("train","validation"),
         xlim=c(5,300),ylim=c(5,300))

### Output cal data for VIP models
write.csv(cal.plsr.data,file=paste(out.dir,'FFT_Leaf_LMA_Calibration_Dataset.csv',
                                   sep=""), row.names=FALSE)

### Standardized model
dims <- dim(cal.plsr.data)
k <- round(dims[1]/10)
segs = cvsegments(dims[1],k = k, type="random")
LeafLMA.pls.stand = plsr(LMA_g_DW_m2~Spectra,scale=TRUE,ncomp=15,validation="CV",segments=segs,
                       trace=TRUE, method = "oscorespls", data=cal.plsr.data)
rm(segs)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# External validation - Use to test initial and final PLSR models
dims <- dim(val.plsr.data)
row.names(val.plsr.data) <- seq(len=nrow(val.plsr.data))
val.spec <- as.matrix(droplevels(val.plsr.data[,17:dims[2]]))
val.plsr.data2 <- data.frame(val.plsr.data[,1:16],Spectra=I(val.spec))
val.plsr.data <- val.plsr.data2
rm(val.plsr.data2)

RMSEP(LeafLMA.pls, newdata = val.plsr.data)
plot(RMSEP(LeafLMA.pls,estimate=c("test"),newdata = val.plsr.data), main="MODEL RMSEP",
     xlab="NUM OF COMPONENTS",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)

R2(LeafLMA.pls, newdata = val.plsr.data)
plot(R2(LeafLMA.pls,estimate=c("test"),newdata = val.plsr.data), main="MODEL R2",
     xlab="NUM OF COMPONENTS",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)

# Quick validation diagnostic plots
predplot(LeafLMA.pls, ncomp = 9:11, newdata = val.plsr.data, asp = 1, line = TRUE,
         which = c("train","validation", "test"),
         xlim=c(5,300),ylim=c(5,300))

### Output validation data for VIP model
write.csv(val.plsr.data,file=paste(out.dir,'FFT_Leaf_LMA_Validation_Dataset.csv',
                                   sep=""), row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Jackknife test here.  Determine optimal number of components
dims = dim(cal.plsr.data)
iterations <- 30
k <- 20
comps <- 15
jk.out <- array(data=NA,dim=c(iterations,comps+1)) # num of components plus intercept
for (i in 1:iterations) {
  print(paste("Iteration: ",i,sep=""))
  segs = cvsegments(dims[1],k = k, type="random")
  LeafLMA.pls = plsr(log10(LMA_g_DW_m2)~Spectra,scale=FALSE,ncomp=comps,validation="CV",segments=segs,
                     trace=TRUE, method = "oscorespls", data=cal.plsr.data)
  vec <- as.vector(RMSEP(LeafLMA.pls)$val)
  vec.sub.adj <- vec[seq(2,length(vec),2)]
  jk.out[i,] <- vec.sub.adj
}

# Boxplot of results
boxplot(jk.out, xaxt="n",xlab="NUM OF COMPONENTS",ylab="RMSEP") 
numcomps <- comps+1
axis(1,at=1:numcomps,0:15)
box(lwd=2.2)

ttest <- t.test(jk.out[,12],jk.out[,13],alternative = "two.sided",paired=F,var.equal = FALSE)
ttest

rm(LeafLMA.pls,comps,i,iterations,vec,vec.sub.adj,ttest,segs,numcomps)
dims <- dim(cal.plsr.data)
k <- round(dims[1]/10)
segs = cvsegments(dims[1],k = k, type="random")
LeafLMA.pls = plsr(LMA_g_DW_m2~Spectra,scale=FALSE,ncomp=15,validation="CV",segments=segs,
                   trace=TRUE, method = "oscorespls", data=cal.plsr.data)

# Examine raw PLSR output
summary(LeafLMA.pls)
plot(RMSEP(LeafLMA.pls), legendpos = "topright")
names(LeafLMA.pls)

predplot(LeafLMA.pls, ncomp = 9:11, asp = 1, line = TRUE,which = c("train","validation"),
         xlim=c(5,300),ylim=c(5,300))
rm(segs)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Calculate Q2 statistic
dims <- dim(cal.plsr.data)
PRESS = LeafLMA.pls$validation$PRESS
SS = sum((cal.plsr.data$LMA_g_DW_m2)^2)
TSS = sum((cal.plsr.data$LMA_g_DW_m2-mean(cal.plsr.data$LMA_g_DW_m2))^2)
Q2 = 1-(PRESS/TSS)  # using SS not TSS

# Calculate RMSECV
RMSECV = PRESS/dims[1]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Generate final model diagnostics
minpress <- which.min(as.vector(PRESS))
minpress
ncomp <- 11 # Final number of components

# Obs validation data
LMA.val <- val.plsr.data$LMA_g_DW_m2
# Predicted val data
pred.val.data <- as.vector(predict(LeafLMA.pls, newdata = val.plsr.data, ncomp=ncomp, 
                                   type="response")[,,1])

# Component selection figure
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Model_Component_Diagnostics.pdf',sep=""),height=8,width=8)

boxplot(jk.out, xaxt="n",xlab="NUM OF COMPONENTS",ylab="RMSEP") 
axis(1,at=1:16,0:15)
box(lwd=2.2)

par(mfrow=c(2,2))
# PRESS
plot(as.vector(PRESS),main="MODEL PRESS",xlab="NUM OF COMPONENTS",ylab="PRESS",cex=1.5,lty=1)
abline(v=ncomp,lty=2,col="black",lwd=2)
abline(v=minpress,lty=2,col="dark grey",lwd=2)
legend("bottomleft",legend=c("Best","Min/Max"),lty=2,col=c("black","dark grey"),lwd=2,,bty = "n")
# RMSEP/RMSECV
plot(RMSEP(LeafLMA.pls,estimate=c("train","CV","test"),newdata = val.plsr.data), main="MODEL RMSEP",
     xlab="NUM OF COMPONENTS",ylab="MODEL RMSEP",lty=c(1,2,2),col=c("black","red","blue"),cex=1.5,lwd=2)
legend("bottomleft",legend=c("Train","CV","Test/Val"),
       col=c("black","red","blue"),lty=c(1,2,2),lwd=c(2,2,2),bty = "n")
abline(v=ncomp,lty=2,col="black",lwd=2)
abline(v=minpress,lty=2,col="dark grey",lwd=2)
# Q2
plot(as.vector(Q2),main="MODEL Q2",xlab="NUM OF COMPONENTS",ylab="Q2",cex=1.5)
abline(v=ncomp,lty=2,col="black",lwd=2)
abline(v=minpress,lty=2,col="dark grey",lwd=2)
# R2 plot
plot(R2(LeafLMA.pls,estimate=c("train","CV","test"),newdata = val.plsr.data), main="MODEL R2",
     xlab="NUM OF COMPONENTS",ylab="MODEL R2",lty=c(1,2,2),col=c("black","red","blue"),cex=1.5,lwd=2)
legend("bottomright",legend=c("Train","CV","Test/Val"),
       col=c("black","red","blue"),lty=c(1,2,2),lwd=c(2,2,2),bty = "n")
abline(v=ncomp,lty=2,col="black",lwd=2)
abline(v=minpress,lty=2,col="dark grey",lwd=2)

# Train/Val/Test component selection plots
predplot(LeafLMA.pls, ncomp = 1:15, newdata = val.plsr.data, asp = 1, line = TRUE,
         which = c("train","validation", "test"),
         xlim=c(5,300),ylim=c(5,300))
dev.off()

# Build output dataset
dims = dim(cal.plsr.data)
n=dims[1]
LMA = cal.plsr.data$LMA_g_DW_m2
cal.plsr.pred = as.vector(LeafLMA.pls$fitted.values[,,ncomp]) # Model fitted values. Predicted values
cal.plsr.CVpred <- as.vector(LeafLMA.pls$validation$pred[,,ncomp]) # CV pred values
cal.residuals = LMA - cal.plsr.pred # Obs minus pred
cal.CVresiduals <- as.vector(LeafLMA.pls$residuals[,,ncomp]) # CV pred residuals
cal.output = data.frame(cal.plsr.data[,1:16],cal.plsr.pred,cal.plsr.CVpred,cal.residuals,cal.CVresiduals)
names(cal.output) = c("Sample_Name","Sample_Year","Site","Plot","Site_Plot","Species",
                      "Height","Age","BL_CON","Water_Perc","LDMC_g_g","EWT_g_cm2","SLA_cm2_g_DW",
                      "SLA_m2_kg_DW","LMA_g_DW_m2","LMA_g_DW_cm2","PLSR_Pred_LMA_g_DW_m2",
                      "PLSR_CV_Pred_LMA_g_DW_m2",
                      "PLSR_LMA_Residuals","PLSR_CV_LMA_Residuals")

# Model statistics
MSECV = mean(cal.CVresiduals^2)
RMSECV = sqrt(MSECV)
PERC_RMSE = (RMSECV/(max(LMA)-min(LMA)))*100
Rsq <- R2(LeafLMA.pls)$val[,,ncomp+1] # Cal/Training R2. Have to add 1, intercept included, i.e. cumulative
Rsq.val <- R2(LeafLMA.pls,newdata = val.plsr.data)$val[,,ncomp+1] # Val R2. Have to add 1, intercept included, i.e. cumulative
Model.bias = mean(cal.plsr.CVpred)-mean(LMA)
names(Model.bias)="Model_bias"
Rsq
Rsq.val
RMSECV
PERC_RMSE
Model.bias

# PLSR Summary statistics
cal_sum_stats = data.frame(Train_Rsq=Rsq,Val_Rsq=Rsq.val,RMSECV=RMSECV,PERC_RMSE=PERC_RMSE,
                           Model_bias=Model.bias)
cal_sum_stats


## PLSR Observed versus predicted plot & independent val plot using withheld samples
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Calibration_Plot.pdf',sep=""),height=12,width=10)
par(mfrow=c(3,2),mar=c(5,5,1,1)) #B, L, T, R
# Cal plot
plot(cal.plsr.pred,LMA,xlim=c(5,250),ylim=c(5,250),pch=21,bg="grey60",
     cex=1.5,xlab="PREDICTED",ylab="OBSERVED",main=paste("Leaf LMA Calibration -- n: ",n))
points(cal.plsr.CVpred,LMA,pch=21,cex=1.5,bg="black")
legend("topleft",legend=c(paste("Train R2 = ",round(Rsq,2)),paste("RMSECV = ",round(RMSECV,2)),
                          paste("Perc RMSECV = ",round(PERC_RMSE,2)),paste("Bias = ",round(Model.bias,4))),
       cex=1.6,bty="n")
abline(0,1,lty=2)
box(lwd=2)

# Val plot
n <- length(LMA.val)
val.residuals <- LMA.val- pred.val.data
MSE.val <- mean(val.residuals^2)
RMSE.val <- sqrt(MSE.val)
Val.bias <- mean(pred.val.data)-mean(LMA.val)
plot(pred.val.data,LMA.val,xlim=c(5,250),ylim=c(5,250),pch=21,bg="grey60",
     cex=1.5,xlab="PREDICTED",ylab="OBSERVED",main=paste("Leaf LMA Validation -- n: ",n))
legend("topleft",legend=c(paste("Val R2 = ",round(Rsq.val,2)),paste("RMSE = ",round(RMSE.val,2)),
                          paste("Bias = ",round(Val.bias,4))),cex=1.6,bty="n")
abline(0,1,lty=2)
box(lwd=2)

# Cal Residuals
plot(LMA,cal.residuals,xlab="LEAF LMA (gDW / m2)",ylab="PLSR Residuals",pch=21,
     bg="grey60",cex=1.5)
abline(h=0,lty=2,col="grey60")
box(lwd=2)

# Val Residuals
plot(LMA.val,val.residuals,xlab="LEAF LMA (gDW / m2)",ylab="PLSR Residuals",pch=21,
     bg="grey60",cex=1.5)
abline(h=0,lty=2,col="grey60")
box(lwd=2)

# Boxplots
hist(cal.residuals)
hist(val.residuals)

dev.off() # End of fig


### Output val dataset
val.output = data.frame(val.plsr.data[,1:16],pred.val.data,val.residuals)
names(val.output) = c("Sample_Name","Sample_Year","Site","Plot","Site_Plot","Species",
                      "Height","Age","BL_CON","Water_Perc","LDMC_g_g","EWT_g_cm2","SLA_cm2_g_DW",
                      "SLA_m2_kg_DW","LMA_g_DW_m2","LMA_g_DW_cm2","PLSR_Pred_LMA_g_DW_m2",
                      "PLSR_LMA_Residuals")

### Scores plot
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Scores_Plot.pdf',sep=""),height=8,width=8)
plot(LeafLMA.pls, plottype = "scores", comps = 1:ncomp)
dev.off()

### Loadings
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Loadings_Plot.pdf',sep=""),height=12,width=10)
par(mfrow=c(2,1))
plot(LeafLMA.pls, plottype = "loadings", comps = 1:4,
     legendpos = "topleft",xlab = "INDEX (500-2400nm)")
plot(LeafLMA.pls, plottype = "loadings", comps = 5:ncomp,
     legendpos = "topleft",xlab = "INDEX (500-2400nm)")
dev.off()

### Loading weights and coefficients
waves = seq(500,2400,1)
weights = loading.weights(LeafLMA.pls)[,1]
coefs = coef(LeafLMA.pls,ncomp=ncomp,intercept=FALSE)
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Loading_Weights_Coeff_Plot.pdf',sep=""),height=12,width=10)
par(mfrow=c(2,1))
plot(weights, lwd=3,xlab = "INDEX (500-2400nm)",cex=0.01)
lines(weights,lwd=3)
abline(h=0,lty=2,lwd=1.5,col="grey60")
plot(waves,coefs,lwd=3,xlab = "WAVELENGTH (nm)", cex=0.01)
lines(waves,coefs,lwd=3)
abline(h=0,lty=2,lwd=1.5,col="grey60")
dev.off()

# VIP Code
source('/Users/serbin/Data/Dropbox/MANUSCRIPTS/FERST_LAB/Dissertation_Manuscripts/R_scripts/VIP.R')

# VIP Plot
waves=seq(500,2400,1)
coefs = coef(LeafLMA.pls,ncomp=ncomp,intercept=FALSE) # WITHOUT INTERCEPT FOR PLOTTING
vips = VIP(LeafLMA.pls)[ncomp,]

pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Coeff_VIP_Plot.pdf',sep=""),height=10,width=10)
par(mfrow=c(2,1))
plot(waves,coefs,cex=0.01,xlab="WAVELENGTH (nm)",ylab="REG COEF")
lines(waves,coefs,lwd=2.5)
abline(h=0,lty=2,col="dark grey")

plot(waves,vips,xlab="WAVELENGTH (nm)",ylab="VIP",cex=0.01)
lines(waves,vips,lwd=3)
abline(h=0.8,lty=2,col="dark grey")
dev.off()
#--------------------------------------------------------------------------------------------------#


#---------------- Export Model Output -------------------------------------------------------------#

# Observed versus predicted
write.csv(cal.output,file=paste(out.dir,'FFT_Leaf_LMA_Observed_PLSR_CV_Pred_',ncomp,
                                'comp.csv',sep=""),row.names=FALSE)

# Validation data
write.csv(val.output,file=paste(out.dir,'FFT_Leaf_LMA_Observed_Val_PLSR_Pred_',ncomp,
                                'comp.csv',sep=""),row.names=FALSE)

# Model coefficients
coefs = coef(LeafLMA.pls,ncomp=ncomp,intercept=TRUE)
write.csv(coefs,file=paste(out.dir,'FFT_Leaf_LMA_PLSR_Coefficients_',ncomp,'comp.csv',sep=""),
          row.names=TRUE)

# standardized
coefs = coef(LeafLMA.pls.stand,ncomp=ncomp,intercept=TRUE)
write.csv(coefs,file=paste(out.dir,'FFT_Leaf_LMA_Standardized_PLSR_Coefficients_',ncomp,'comp.csv',sep=""),
          row.names=TRUE)

# Model loading weights
write.csv(weights,file=paste(out.dir,'FFT_Leaf_LMA_PLSR_Loading_Weights_Comp1.csv',sep=""))

# PLSR VIP
write.csv(vips,file=paste(out.dir,'FFT_Leaf_LMA_PLSR_VIPs_',ncomp,'comp.csv',sep=""))

# PLSR Model stats
write.csv(cal_sum_stats,file=paste(out.dir,'FFT_Leaf_LMA_PLSR_Statistics_',ncomp,'comp.csv',
                                   sep=""), row.names=FALSE)

# Remove temp objects before proceeding
rm(jk.out,Q2,PRESS,cal_sum_stats,LMA,LMA.val,MSE.val,MSECV,Model.bias,
   PERC_RMSE,RMSE.val,RMSECV,Rsq,Rsq.val,SS,TSS,Val.bias,cal.CVresiduals,cal.plsr.CVpred,
   cal.plsr.pred,cal.residuals,coefs,dims,minpress,n,pred.val.data,val.residuals,vips,
   waves,weights)
#--------------------------------------------------------------------------------------------------#


#---------------- Jackknife model validation ------------------------------------------------------#
ncomp
ncomp <- ncomp  # Determined previously
resamples <- 1000 #1000
output.jackknife.stats <- data.frame(Rsq=rep(NA,resamples),RMSEP=rep(NA,resamples),
                                     PERC_RMSEP=rep(NA,resamples), Bias=rep(NA,resamples))
output.jackknife.coefs <- array(data=NA,dim=c(resamples,
                               dim(coef(LeafLMA.pls,ncomp=ncomp,intercept=TRUE))[1]))

for (i in 1:resamples) {
  rows = sample(1:nrow(cal.plsr.data),floor(0.7*nrow(cal.plsr.data)))
  cal.data = cal.plsr.data[rows,]
  val.data = cal.plsr.data[-rows,]
  
  dimsCal = dim(cal.data)
  dimsVal = dim(val.data)
  
  ### Build PLSR model with training data
  LeafLMA.pls.jack <- plsr(LMA_g_DW_m2~Spectra, ncomp=ncomp, SCALE=FALSE, validation="none", data=cal.data)
  
  ### Estimate independent (Validation) samples
  LMA.val <- val.data$LMA_g_DW_m2
  pred.val.data = as.vector(predict(LeafLMA.pls.jack,newdata=val.data$Spectra,
                                    ncomp=ncomp,type="response")[,,1])
  output.jackknife.coefs[i,] <- as.vector(coef(LeafLMA.pls.jack,ncomp=ncomp,intercept=TRUE))
  
  # Error statistics
  n <- length(LMA.val)
  Rsq.val <- R2(LeafLMA.pls.jack,newdata = val.data)$val[,,ncomp+1] 
  RMSEP.val <- RMSEP(LeafLMA.pls.jack,newdata = val.data)$val[,,ncomp+1]
  PERC_RMSEP = (RMSEP.val/(max(LMA.val)-min(LMA.val)))*100
  val.residuals <- LMA.val- pred.val.data
  #  MSE.val <- mean(val.residuals^2)
  #  RMSE.val <- sqrt(MSE.val)
  Val.bias <- mean(pred.val.data)-mean(LMA.val)
  
  ### Store results of iteration i
  output.jackknife.stats[i,1] = Rsq.val
  output.jackknife.stats[i,2] = RMSEP.val
  output.jackknife.stats[i,3] = PERC_RMSEP
  output.jackknife.stats[i,4] = Val.bias
  
  print(paste("Running Iteration",i))
  print(paste("Stats: ","Rsq ",round(Rsq.val,2)," / RMSEP ",round(RMSEP.val,2), " / %RMSEP ",
              round(PERC_RMSEP,2)," / Bias ",round(Val.bias,2),sep="" ) )
  flush.console()  # force the output
  
  # Remove temp objects
  rm(cal.data,val.data,n,LeafLMA.pls.jack,LMA.val,Rsq.val,pred.val.data,RMSEP.val,PERC_RMSEP,
     val.residuals,Val.bias)
}
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Generate jackknife stats object for output
LeafLMA_medRSQ <- median(output.jackknife.stats[,1])
LeafLMA_meanRSQ <- mean(output.jackknife.stats[,1])
LeafLMA_sdRSQ <- sd(output.jackknife.stats[,1])
LeafLMA_minRSQ <- min(output.jackknife.stats[,1])
LeafLMA_maxRSQ <- max(output.jackknife.stats[,1])

LeafLMA_medRMSE <- median(output.jackknife.stats[,2])
LeafLMA_meanRMSE <- mean(output.jackknife.stats[,2])
LeafLMA_sdRMSE <- sd(output.jackknife.stats[,2])
LeafLMA_minRMSE <- min(output.jackknife.stats[,2])
LeafLMA_maxRMSE <- max(output.jackknife.stats[,2])

LeafLMA_medBias <- median(output.jackknife.stats[,4])
LeafLMA_meanBias <- mean(output.jackknife.stats[,4])
LeafLMA_sdBias <- sd(output.jackknife.stats[,4])
LeafLMA_minBias <- min(output.jackknife.stats[,4])
LeafLMA_maxBias <- max(output.jackknife.stats[,4])

summary.statistics <- data.frame(MEDIAN_RSQ=LeafLMA_medRSQ,MEAN_RSQ=LeafLMA_meanRSQ,SD_RSQ=LeafLMA_sdRSQ,
                                 MIN_RSQ=LeafLMA_minRSQ,MAX_RSQ=LeafLMA_maxRSQ,MEDIAN_RMSE=LeafLMA_medRMSE,
                                 MEAN_RMSE=LeafLMA_meanRMSE,SD_RMSE=LeafLMA_sdRMSE,MIN_RMSE=LeafLMA_minRMSE,
                                 MAX_RMSE=LeafLMA_maxRMSE,MEDIAN_BIAS=LeafLMA_medBias,
                                 MEAN_BIAS=LeafLMA_meanBias,SD_BIAS=LeafLMA_sdBias,MIN_BIAS=LeafLMA_minBias,
                                 MAX_BIAS=LeafLMA_maxBias)
summary.statistics
write.csv(summary.statistics,file=paste(out.dir,'FFT_Leaf_LMA_Jackkife_Summary_Stats.csv',sep=""),
          row.names=FALSE)
rm(i,dimsCal,dimsVal,rows,LeafLMA_medRSQ,LeafLMA_meanRSQ,LeafLMA_sdRSQ,LeafLMA_minRSQ,LeafLMA_maxRSQ,
   LeafLMA_medRMSE,LeafLMA_meanRMSE,LeafLMA_sdRMSE,LeafLMA_minRMSE,LeafLMA_maxRMSE,LeafLMA_medBias,
   LeafLMA_meanBias,LeafLMA_sdBias,LeafLMA_minBias,LeafLMA_maxBias,summary.statistics)
#----------------------------------------------------------------------------------------------------#


#---------------- Jackknife diagnostic histogram figure ---------------------------------------------#
# Histogram statistics
rsqHist = hist(output.jackknife.stats[,1],plot=FALSE)
rmseHist = hist(output.jackknife.stats[,2],plot=FALSE)
biasHist = hist(output.jackknife.stats[,4],plot=FALSE)

rsqHist_sum = sum(rsqHist$counts)
rsqHist_perc = (rsqHist$counts/rsqHist_sum)*100 

rmseHist_sum = sum(rmseHist$counts)
rmseHist_perc = (rmseHist$counts/rmseHist_sum)*100 

biasHist_sum = sum(biasHist$counts)
biasHist_perc = (biasHist$counts/biasHist_sum)*100

LeafLMA_RSQ_mids = rsqHist$mids
LeafLMA_RMSE_mids = rmseHist$mids
LeafLMA_BIAS_mids = biasHist$mids

# Generate figure
pdf(file=paste(out.dir,"FFT_Leaf_LMA_Jackknife_PLSR_Histogram.pdf",sep=""),width= 20, height= 7)
par(mfrow=c(1,3),mar=c(6,5,1,0.3),oma=c(0,0,0,0))
barplot(rsqHist_perc, border = "dark grey",xlab =expression(paste("Model ",R^2)),
        ylim=c(0,range(rsqHist_perc)[2]+4),ylab = "Percent (%)",axes = TRUE, axisnames = TRUE,
        axis.lty = 1,lwd=2.5,cex.axis=1.7,cex.names=1.7,names.arg=LeafLMA_RSQ_mids, cex.lab=2.5)
box(lwd=2.5)

# RMSE
barplot(rmseHist_perc, border = "dark grey",ylim=c(0,range(rmseHist_perc)[2]+4),xlab="Model RMSE",
        ylab =,axes = TRUE, axisnames = TRUE,axis.lty = 1,lwd=2.5,cex.axis=1.7,cex.names=1.7,
        names.arg=LeafLMA_RMSE_mids, cex.lab=2.5)
box(lwd=2.5)

# Model bias
barplot(biasHist_perc, border = "dark grey",ylim=c(0,range(biasHist_perc)[2]+4),xlab="Model Bias",
        ylab =,axes = TRUE, axisnames = TRUE,axis.lty = 1,lwd=2.5,cex.axis=1.7,cex.names=1.7,
        names.arg=LeafLMA_BIAS_mids, cex.lab=2.5)
box(lwd=2.5)

dev.off()

rm(rsqHist,rmseHist,biasHist,rsqHist_sum,rsqHist_perc,rmseHist_sum,rmseHist_perc,biasHist_sum,
   biasHist_perc,LeafLMA_RSQ_mids,LeafLMA_RMSE_mids,LeafLMA_BIAS_mids)
#--------------------------------------------------------------------------------------------------#


#---------------- Output jackknife results --------------------------------------------------------#
write.csv(output.jackknife.stats,file=paste(out.dir,'FFT_Leaf_LMA_Jackkife_PLSR_Resutls.csv',sep=""),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# PLSR Coef Jackknife Diagnostics
dims <- dim(output.jackknife.coefs)
boxplot(output.jackknife.coefs[,1])
plot.coefs <- output.jackknife.coefs[,2:dims[2]]
plot.min <- min(output.jackknife.coefs[,2:dims[2]])
plot.max <- max(output.jackknife.coefs[,2:dims[2]])

jk.intercepts <- output.jackknife.coefs[,1]

# Stats
coef.means <- colMeans(plot.coefs)
sd.coef <- apply(plot.coefs,MARGIN=2,FUN=function(x)sd(x))
min.coef <- apply(plot.coefs,MARGIN=2,FUN=function(x)min(x))
max.coef <- apply(plot.coefs,MARGIN=2,FUN=function(x)max(x))
coef.quant <- apply(plot.coefs,2,quantile,probs=c(0.05,0.95))
intercepts.quant <- quantile(jk.intercepts,probs=c(0.05,0.95))

# T Test
x <- resamples # From Jackknife above
#results <- apply(plot.coefs,2,t.test,alternative = c("two.sided"))
results <- apply(plot.coefs,2, function(plot.coefs) {
  t.test(x = plot.coefs[1:x])$p.value}) 


# JK Coef Figure
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Jackknife_Coef_Diagnostic_Figure.pdf',sep=""),height=20,width=10)
par(mfrow=c(3,1),mar=c(4,4,1,1.2)) #B,L,T,R

waves <- seq(500,2400,1)
for (i in 1:resamples){
  if (i==1){
    plot(waves,plot.coefs[i,],type="l",lwd=1.3,xlab="Wavelength (nm)",
         ylab="PLSR Coefficients",ylim=c(plot.min,plot.max))
    abline(h=0,lty=2,col="grey",lwd=1.5)
  } else {
    lines(waves,plot.coefs[i,],lwd=1.3)
  }
}

coefs = as.vector(coef(LeafLMA.pls,ncomp=ncomp,intercept=FALSE))
plot(waves,coefs,type="l",lwd=4,ylim=c(plot.min,plot.max))

# Min/Max
polygon(c(waves ,rev(waves)),c(max.coef, rev(min.coef)),col="grey50",border=NA)
lines(waves,min.coef,lty=1,lwd=3,col="grey50")
lines(waves,max.coef,lty=1,lwd=3,col="grey50")

# 95% CIs
polygon(c(waves ,rev(waves)),c(coef.quant[2,], rev(coef.quant[1,])),col="grey70",border=NA)
lines(waves,coef.quant[1,],lty=1,lwd=2,col="grey70")
lines(waves,coef.quant[2,],lty=1,lwd=2,col="grey70")

# replot the mean and zero line
lines(waves,coefs,lwd=4)
abline(h=0,lty=2,col="grey",lwd=1.5)

legend("bottomright",legend=c("Mean","95% CI"),lty=c(1,1),
       col=c("black","dark grey"),lwd=3)
box(lwd=2.2)

# P-vals
plot(waves,results,pch=21,bg="grey80",ylab="P-value",xlab="Wavelength (nm)",
     cex=4)

# Coeff trace plot
par(mfrow=c(3,1),mar=c(4,4,2,1.2)) #B,L,T,R
for (i in 1:length(waves)){
  wavelength <- waves[i]
  median <- median(plot.coefs[,i])
  mean <- mean(plot.coefs[,i])
  plot(plot.coefs[,i],type="l",main=paste(wavelength),xlab="Iteration",
       ylab="PLSR Coefficients")
  abline(h=median,lty=2,col="grey50",lwd=1.4)
  #abline(h=mean,lty=2,col="grey50",lwd=1.4)
}

dev.off()

#--------------------------------------------------------------------------------------------------#


#---------------- Output jackknife results --------------------------------------------------------#
# JK Coefficents
out.jk.coefs <- data.frame(Iteration=seq(1,resamples,1),jk.intercepts,plot.coefs)
names(out.jk.coefs) <- c("Iteration","Intercept",paste("Wave_",seq(500,2400,1),sep=""))
write.csv(out.jk.coefs,file=paste(out.dir,'FFT_Leaf_LMA_Jackkife_PLSR_Coefficients.csv',sep=""),
          row.names=FALSE)

# Coeff quantiles
out.coef.quant <- array(data=NA,dim=c(2,1903))
out.coef.quant[1,1] <- "5%"
out.coef.quant[2,1] <- "95%"
out.coef.quant[1,2] <- intercepts.quant[[1]]
out.coef.quant[2,2] <- intercepts.quant[[2]]
out.coef.quant[,3:1903] <- coef.quant
out.coef.quant <- data.frame(out.coef.quant)

names(out.coef.quant) <- c("Quantile","Intercept",paste("Wave_",waves,sep=""))

write.csv(out.coef.quant,file=paste(out.dir,'FFT_Leaf_LMA_Jackkife_PLSR_Coefficient_Quantiles.csv',sep=""),
          row.names=TRUE)
# P-vals
out.pvals <- data.frame(Wavelength=paste("Wave_",seq(500,2400,1),sep=""),Pval=results)
write.csv(out.pvals,file=paste(out.dir,'FFT_Leaf_LMA_Jackkife_PLSR_Coefficient_Pvals.csv',sep=""),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# JK Val plot
dims <- dim(output.jackknife.coefs)
intercepts <- output.jackknife.coefs[,1]
jk.coef.test.output <- array(data=NA,dim=c(dim(val.spec)[1],dims[1]))
for (i in 1:length(intercepts)){
  #for (i in 1:10){
  print(paste("Iteration: ",i,sep=""))
  coefs <- as.vector(output.jackknife.coefs[i,2:dims[2]])
  temp <- val.spec %*% coefs  # Updated: Using matrix mult.
  vals = data.frame(rowSums(temp))+intercepts[i]
  jk.coef.test.output[,i] <- vals[,1]
}

# Figure 1
# PLSR output: val.output
# 95% CI data: jk.coef.test.output
pred.quant <- apply(jk.coef.test.output,1,quantile,probs=c(0.05,0.95))
pred.quant.ll <- pred.quant[1,]
pred.quant.ul <- pred.quant[2,]

obs <- as.matrix(val.output$LMA_g_DW_m2)
pred <- as.matrix(val.output$PLSR_Pred_LMA_g_DW_m2)
y <- as.vector(val.output$LMA_g_DW_m2)
x <- as.vector(val.output$PLSR_Pred_LMA_g_DW_m2)

reg <- lm(y~x)
summary(reg)
new <- data.frame(x=seq(-20,400,0.1))
predint <- predict(reg, new, interval="prediction")
confint <- predict(reg, new, interval="confidence")
int <- summary(reg)$coefficients[1]

pred.rsq <- unlist(unlist(R2(LeafLMA.pls, newdata = val.plsr.data))[ncomp+1])
pred.rmsep <- unlist(unlist(RMSEP(LeafLMA.pls, newdata = val.plsr.data))[ncomp+1])
residuals <- pred - obs
P.RMSEP <- sqrt(colMeans(residuals^2)) # should be similar to other

# Obs vs Pred
pdf(paste(out.dir,'FFT_Leaf_LMA_PLSR_Validation_Diagnostic_Figure.pdf',sep=""),height=8.5,width=8.8)
par(mar=c(4.0,4.4,1,0.9)) #B,L,T,R
plotCI(val.output$PLSR_Pred_LMA_g_DW_m2,val.output$LMA_g_DW_m2,li=pred.quant.ll,gap=0.009,sfrac=0.004,lwd=1.6,
       ui=pred.quant.ul,err="x",pch=21,col="black",pt.bg="grey70",xlim=c(5,250),cex=1.3,
       ylim=c(5,250),xlab="Predicted",ylab="Observed",main="Leaf LMA (gDW / m2)",cex.axis=1.5,cex.lab=1.8)
abline(0,1,lty=2,lw=2)
lines(new$x,predint[,1])
lines(new$x,confint[,2],lwd=2,col="grey")
lines(new$x,confint[,3],lwd=2,col="grey")
lines(new$x,predint[,2],lwd=2.5)
lines(new$x,predint[,3],lwd=2.5)
legend("topleft",legend=c(paste("R2 = ",round(pred.rsq,2),sep=""),paste("RMSE = ",round(pred.rmsep,2),sep=""),
                          paste("Reg. Bias =", round(int,2))),bty="n",cex=1.7)
legend("bottomright",legend=c("1:1 Line","Regression Line","95% Conf. Int.","95% Pred. Int."),
       lty=c(2,1,1,1),lwd=c(1.6,1.3,2.1,2.1),col=c("black","black","grey","black"))
box(lwd=2.2)
dev.off()

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Output JK Coefficient test results
jk.coef.test.output2 <- data.frame(jk.coef.test.output)
names(jk.coef.test.output2) <- paste("Iteration.",seq(1,resamples,1),sep="")
jk.coef.test.output2 <- data.frame(Observed.Values=val.output$LMA_g_DW_m2, jk.coef.test.output2)

write.csv(jk.coef.test.output2,file=paste(out.dir,'FFT_Leaf_LMA_Jackkife_PLSR_Val_Data_Output.csv',sep=""),
          row.names=FALSE)

write.csv(predint,file=paste(out.dir,'FFT_Leaf_LMA_PLSR_Val_Prediction_Intervals.csv',sep=""),
          row.names=FALSE)

write.csv(confint,file=paste(out.dir,'FFT_Leaf_LMA_PLSR_Val_Confidence_Intervals.csv',sep=""),
          row.names=FALSE)
#--------------------------------------------------------------------------------------------------#

# End of program