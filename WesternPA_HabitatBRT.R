###################################################################################################
## CUSTOM CODE FOR ANALYZING PENNSYLVANIA WILDS HABITAT STRUCTURE DATA
## AUTHOR: Nicole Michel, National Audubon Society
## DATE: February 2016
## 
###################################################################################################


# set working directory
setwd("~/GitHub/Climate2/western_PA_density")

library(dismo)

# read in source code with functions
source("WesternPA_Density_Functions.R")

# define the path to where your input and output data files are stored
pathtofiles <- "C:/Users/nmichel/Box Sync/2_Projects/Habitats/Forests/Western PA Restoration Strategy/"

# read in data
ptctTot <- read.csv(file=paste(pathtofiles,"WPAData_TotalBySpeciesPointSurvey.csv", sep=""))
ptctTot.spp <- subset(ptctTot, ptctTot$elem_name=="BTBW")
offsets <- read.csv(file=paste(pathtofiles,"DetOffset_BTBW_2015.csv", sep=""))
brtdata <- merge(ptctTot.spp, offsets[,c("pt_id","DetOffset")], by.x="pt_id", by.y="pt_id")

str(brtdata)

# subset to include only the pt_id, TotSing (count), offset, and those predictors you want to include in the model
brtdata <- brtdata[,c("pt_id","TotSing","elev","slope_per","aspect","topo_pos","bedrock","litter_duf",
                          "bryophyte_","lg_rocks","wood","other","sm_rocks","water","sand","bare_soil","leaf_phen","leaf_type","phys_cl",
                          "pa_comm_code","basal_area","spec_t1_co", "spec_t1_ht", "spec_t2_co", "spec_t2_ht", 
                          "spec_t3_co", "spec_t3_ht", "spec_s1_co", "spec_s1_ht", "spec_s2_co",  "spec_s2_ht" , 
                          "spec_h_cov", "spec_h_htc", "spec_n_cov", "spec_n_htc", "spec_v_cov", "spec_v_htc", "snags", 
                          "cavity", "water_pres" , "invasive_s", "invasive_1", "rd_pav", "rd_unpaved", "utility_li", "trail_pave", "trail_grav", 
                          "trail_unim", "ditch_lg", "ditch_sm", "grading", "equip_trac", "pine_plant", "clearcut", "logging", "mowing", "grazing", 
                          "under_remo" , "herbivory", "fire" , "blowdown", "tree_disea", "pest", "landslide" , "garbage", "gas_devel_", "dist_index","DetOffset","LongFold"),]
# do not include species variables - brts tend to allocate too much weight to factors with lots of levels
# "spec_t1","spec_t2", "spec_t3" , "spec_s1", "spec_s2" , "spec_h", "spec_n", "spec_v", 

##################################################################################
### EXPLORE AND PREP DATA
##################################################################################

# investigate collinearity among continuous predictors. Typical cutoff for GLMs: 0.7. With BRTs I go up to ~0.8-0.9 before worrying

# function to extract only continuous predictors
subset_colclasses <- function(DF, colclasses="numeric") {
  DF[,sapply(DF, function(vec, test) class(vec) %in% test, test=colclasses)]
}
brtcontpreds <- subset_colclasses(brtdata[,c(3:67)], c("numeric", "integer"))
cormat <- cor(brtcontpreds)
max(cormat[which(!is.na(cormat) & cormat[,]<1)])
min(cormat[which(!is.na(cormat))])
cormat


##################################################################################
### BEGIN BRT ANALYSIS
##################################################################################

#### @@@@ USER INPUT NEEDED HERE @@@@ - check to make sure there are no blank (X, X.1, X.2...) columns at the beginning
### TotSing should be in second column. Otherwise model will explode
head(brtdata)
#brtdata$X <- NULL

# add a random number column to the data for later model simplification (will remove predictors that explain no more variation than a random number)
brtdata$RandNum <- sample(c(1:1000), size=nrow(brtdata), replace=T)

# start with lr=0.1, tc=1, bagf=0.5
wpabrt.lr001.bf5.tc1 <- gbm.step(data=brtdata, gbm.x = 3:76, gbm.y = 2, offset=brtdata$DetOffset,
                            family="poisson", fold.vector=brtdata$LongFold,  tree.complexity = 1,
                            learning.rate = 0.1, bag.fraction = 0.5)

# remove columns with no variation or only missing values, then rerun
brtdata$spec_t1_ht <- NULL
brtdata$spec_v_cov <- NULL
brtdata$rd_pav <- NULL
brtdata$rd_unpaved <- NULL
brtdata$trail_pave <- NULL
brtdata$trail_grav <- NULL
brtdata$utility_li <- NULL
brtdata$ditch_lg <- NULL
brtdata$ditch_sm <- NULL
brtdata$grading <- NULL
brtdata$equip_trac <- NULL
brtdata$pine_plant <- NULL
brtdata$mowing <- NULL
brtdata$clearcut <- NULL
brtdata$fire <- NULL
brtdata$other <- NULL
brtdata$pest <- NULL
brtdata$landslide <- NULL
brtdata$garbage <- NULL
brtdata$gas_devel_ <- NULL
brtdata$sand <- NULL

# fit the first model with arbitrarily chosen parameters. 
# Tune learning rate, continuing to lower it until the model fits >1000 trees.
wpabrt.lr003.bf5.tc1 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset,
                               family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 1,
                               learning.rate = 0.003, bag.fraction = 0.5)
# 1650 trees, mean total deviance explained = 1-(0.762/0.937) = 18.7%, cv corr= .31 ## summarize model fit statistics here
# just ignore warnings about offsets, they don't influence model fit

## LEARNING RATE DETERMINED, NOW DETERMINE BAG FRACTION. USUALLY BETWEEN 0.5 - 0.75
##  Arbitrarily increase the bf by 0.05-0.1 and compare with the previous model based on total deviance explained
##   and cv corr. Select the bf that produces the highest total deviance explained and/or cv corr 
##   (if the two measures disagree, go with deviance explained)

wpabrt.lr003.bf55.tc1 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset,
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 1,
                                 learning.rate = 0.003, bag.fraction = 0.55)
# 1550 trees, mean total deviance explained = 1-(0.756/0.937) = 19.3%, cv corr= .313 

wpabrt.lr003.bf6.tc1 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset,
                                  family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 1,
                                  learning.rate = 0.003, bag.fraction = 0.6)
# 1450 trees, mean total deviance explained = 1-(0.754/0.937) = 19.5%, cv corr= .312 = BEST

wpabrt.lr003.bf65.tc1 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset,
                                  family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 1,
                                  learning.rate = 0.003, bag.fraction = 0.65)
# 1750 trees, mean total deviance explained = 1-(0.76/0.937) = 18.9%, cv corr= .312 



## NOW TRY INCREASING TREE COMPLEXITY OF BEST MODEL
wpabrt.lr001.bf6.tc2 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                   family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 2,
                                   learning.rate = 0.001, bag.fraction = 0.6)
# 3350 trees, mean total deviance explained = 1-(0.736/0.937) = 21.5%, cv corr= .306
# we have to reduce the learning rate to 0.001 to get >1000 trees, so need to refit model with lr=0.001 and tc=1 for comparison

wpabrt.lr001.bf6.tc1 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 1,
                                 learning.rate = 0.001, bag.fraction = 0.6)
# 4300 trees, mean total deviance explained = 1-(0.772/0.937) = 17.6%, cv corr= .314

wpabrt.lr001.bf6.tc3 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 3,
                                 learning.rate = 0.001, bag.fraction = 0.6)
# 2350 trees, mean total deviance explained = 1-(0.732/0.937) = 21.9%, cv corr= .299

wpabrt.lr001.bf6.tc4 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 4,
                                 learning.rate = 0.001, bag.fraction = 0.6)
# 1800 trees, mean total deviance explained = 1-(0.731/0.937) = 22.0%, cv corr= .297

wpabrt.lr001.bf6.tc5 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 5,
                                 learning.rate = 0.001, bag.fraction = 0.6)
# 2000 trees, mean total deviance explained = 1-(0.697/0.937) = 25.6%, cv corr= .296

wpabrt.lr001.bf6.tc6 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 6,
                                 learning.rate = 0.001, bag.fraction = 0.6)
# 1950 trees, mean total deviance explained = 1-(0.68/0.937) = 27.4%, cv corr= .29 = BEST

wpabrt.lr001.bf6.tc7 <- gbm.step(data=brtdata, gbm.x = c(3:46,49), gbm.y = 2, offset=brtdata$DetOffset, 
                                 family="poisson", fold.vector=brtdata$LongFold, tree.complexity = 7,
                                 learning.rate = 0.001, bag.fraction = 0.6)
# 1950 trees, mean total deviance explained = 1-(0.68/0.937) = 27.4%, cv corr= .287

spec.best <- wpabrt.lr001.bf6.tc7
save(spec.best, file=paste(pathtofiles,"BTBW_FullBRTModel_2015.rda",sep=""))


#################### SIMPLIFY AND REVIEW BEST MODEL ###################################
summary(spec.best)

# get variables with higher VI than Random
RandomVI <- spec.best$contributions$rel.inf[which(rownames(spec.best$contributions)=="RandNum")]
specSigVars <- subset(spec.best$contributions, spec.best$contributions$rel.inf>RandomVI)
spec.varnames <- rownames(specSigVars)
spec.varno <- which(colnames(brtdata) %in% spec.varnames)

# rerun best BRT model with only those predictors that made it through the simplify process 
spec.simp <- gbm.step(data=brtdata, gbm.x = spec.varno, gbm.y = which(colnames(brtdata)=="TotSing"), 
                      offset=brtdata$DetOffset, family="poisson", fold.vector=brtdata$LongFold, 
                      tree.complexity = spec.best$interaction.depth, learning.rate = 0.001, bag.fraction = spec.best$bag.fraction)
# 1400 trees, mean total deviance explained = 1-(0.721/0.937) = 23.1%, cv corr= .284


## THIS IS YOUR FINAL MODEL. THE STATISTICS TO REPORT ARE: 
DevExpl <- 1-(spec.simp$self.statistics$mean.resid/spec.simp$self.statistics$mean.null)  # Deviance explained
CVcorr <- spec.simp$cv.statistics$correlation.mean  # correlation between observations and predictions. Cross-validated with spatially stratified testing and evaluation datasets

# save simplified model
save(spec.simp, file=paste(pathtofiles,"BTBW_SimplifiedBRTModel_2015.rda", sep=""))

# look at variable importance in final, simplified model
summary(spec.simp)

# produce marginal effects plots for top 3 predictors (just for preliminary viewing, will produce nicer plots below)
gbm.plot22(spec.simp, variable.no=which(spec.simp$var.names=="pa_comm_code"), title.text = "BTBW: PA Community Type", x.label = "PA Community Type", y.label = "Marginal effect on BTBW density", write.title=FALSE, common.scale=TRUE, rug=FALSE, show.contrib=FALSE, smooth=FALSE)
gbm.plot22(spec.simp, variable.no=which(spec.simp$var.names=="aspect"), title.text = "BTBW: Aspect", x.label = "Aspect", y.label = "Marginal effect on BTBW density", write.title=FALSE, common.scale=TRUE, rug=FALSE, show.contrib=FALSE, smooth=FALSE)
gbm.plot22(spec.simp, variable.no=which(spec.simp$var.names=="elev"), title.text = "BTBW: Elevation", x.label = "Elevation", y.label = "Marginal effect on BTBW density", write.title=FALSE, common.scale=TRUE, rug=FALSE, show.contrib=FALSE, smooth=FALSE)
gbm.plot22(spec.simp, variable.no=which(spec.simp$var.names=="spec_s1_co"), title.text = "BTBW: Spec S1 co", x.label = "Spec S1 co", y.label = "Marginal effect on BTBW density", write.title=FALSE, common.scale=TRUE, rug=FALSE, show.contrib=FALSE, smooth=FALSE)
gbm.plot22(spec.simp, variable.no=which(spec.simp$var.names=="slope_per"), title.text = "BTBW: Slope", x.label = "Slope", y.label = "Marginal effect on BTBW density", write.title=FALSE, common.scale=TRUE, rug=FALSE, show.contrib=FALSE, smooth=FALSE)
gbm.plot22(spec.simp, variable.no=which(spec.simp$var.names=="spec_s2_co"), title.text = "BTBW: Spec S2 co", x.label = "Spec S2 co", y.label = "Marginal effect on BTBW density", write.title=FALSE, common.scale=TRUE, rug=FALSE, show.contrib=FALSE, smooth=FALSE)



########################## FIND AND PLOT INTERACTIONS ################

# find interactions
spec.int <- gbm.interactions(spec.simp)
save(spec.int, file=paste(pathtofiles, "BTBW_Interactions.rda", sep=""))

# load(file=paste(pathtofiles, "BTBW_Interactions.rda", sep="")) # this code loads the model output back in

# produce ranked list of most important interactions
head(spec.int$rank.list)

# plot pairwise interactions. Need to update z value so it's just above the maximum value (reported when you run the command)
par(mfrow=c(1,1))
gbm.perspec(spec.simp, which(spec.simp$var.names=="pa_comm_code"), which(spec.simp$var.names=="elev"), leg.coords=c(1,1), z.range=c(0,0.7), x.label="PA Community Type", y.label="Elevation", z.label="Marginal effect on BTBW density") #, z.range=c(0,0.6))


### PLOT VARIABLE IMPORTANCE
# resort VI in increasing order
specVI <- data.frame(spec.simp$contributions)
specVIsort <- specVI[order(specVI$rel.inf),]

plot.new()
jpeg(file=paste(pathtofiles, "BTBW_VariableImportance.jpeg", sep=""), width = 4, height = 4, units = 'in', res = 600)
par((mfrow=c(1,1)), mar=c(5,7,9,2) )
barplot(specVIsort$rel.inf, width=1, horiz=TRUE, beside=TRUE, density=NA, space=1.25, xlim=c(0,(round(max(specVIsort$rel.inf)*1.1))), ylim=c(0,(nrow(specVIsort))),
        names.arg=specVIsort$var, las=1, cex.names = .75, cex.axis = .75, col="blue", border=NA, xlab="Variable Importance") # las=1 rotates lables 90d
dev.off()


############## CHECK TO SEE IF THERE IS SUBSTANTIAL SPATIAL AUTOCORRELATION IN THE MODEL RESIDUALS

# merge in lat/long data (X/Y)
ptctTot <- read.csv(file=paste(pathtofiles,"WPAData_TotalBySpeciesPointSurvey.csv", sep=""))
ptctTot.spp <- subset(ptctTot, ptctTot$elem_name=="BTBW")
brtdataXY <- merge(brtdata, ptctTot.spp [,c("pt_id","X","Y")], by.x="pt_id", by.y="pt_id")

# associate residuals with spatial location. Do only for 2011 (can't have co-located data)
gbm.resid<-residuals(spec.simp,type="pearson")
gbmresxy<-data.frame(gbm.resid,brtdataXY$X,brtdataXY$Y)
colnames(gbmresxy) <- c("Resid","X","Y")

# generate distance matrix, invert it, and replace diagonals with 0
gbmres.dists <- as.matrix(dist(cbind(gbmresxy$X, gbmresxy$Y)))
gbmres.dists.inv <- 1/gbmres.dists
diag(gbmres.dists.inv) <- 0
gbmres.dists.inv[is.infinite(gbmres.dists.inv)] <- 0

# calculate Moran's I of residuals. I = observed. I > 0.3 is evidence of problematic spatial autocorrelation (Lichstein et al., 2002)
library(ape)
Moran.I(gbmresxy$Resid, gbmres.dists.inv) # i = 0.10, p = 0

## plot spatial distribution of residuals. You want the bubbles of different sizes and colors to be ~ randomly distributed
coordinates(gbmresxy)<-c("X","Y")
bubble(gbmresxy,"Resid",col=c("black","darkgrey"),pch=16,maxsize=2,main="BTBW model residuals",identify=FALSE,xlab="Longitude",ylab="Latitude",scales=list(draw=TRUE))


################### PRODUCE PLOTS SHOWING RESPONSE OF PREDICTED DETECTION-CORRECTED DENSITY TO PREDICTORS ###############
###################  (USE INSTEAD OF MARGINAL RESPONSE PLOTS, ABOVE)

######################################################################
## BOOTSTRAP STANDARD ERRORS

# create new data frames for each predictor in the simplified model, with 
#   representative values of the focal predictor and all others held to their means.
#   use mean offset for each level of the focal predictor in the data frame
## PA community code
paccoffset <- vector()
for (i in 1:length(unique(brtdata$pa_comm_code))){
  paccoffset[i] <- mean(brtdata$DetOffset[which(brtdata$pa_comm_code==unique(brtdata$pa_comm_code)[i])])
  if (is.nan(paccoffset[i])){    paccoffset[i] <- mean(brtdata$DetOffset)   }
}
newdata.pacc <- data.frame(pa_comm_code=unique(brtdata$pa_comm_code), 
                             aspect=rep("NE",length(unique(brtdata$pa_comm_code))), 
                             elev=rep(mean(brtdata$elev),length(unique(brtdata$pa_comm_code))), 
                             spec_s1_co=rep(mean(brtdata$spec_s1_co),length(unique(brtdata$pa_comm_code))),
                             slope_per=rep(mean(brtdata$slope_per),length(unique(brtdata$pa_comm_code))),
                             spec_s2_co=rep(mean(brtdata$spec_s2_co),length(unique(brtdata$pa_comm_code))), 
                             DetOffset=paccoffset)
## aspect
aspoffset <- vector()
for (i in 1:length(unique(brtdata$aspect))){
  aspoffset[i] <- mean(brtdata$DetOffset[which(brtdata$aspect==unique(brtdata$aspect)[i])])
  if (is.nan(aspoffset[i])){    aspoffset[i] <- mean(brtdata$DetOffset)   }
}
newdata.aspect <- data.frame(pa_comm_code=rep("CS",length(unique(brtdata$aspect))), 
                             aspect=unique(brtdata$aspect),  
                           elev=rep(mean(brtdata$elev),length(unique(brtdata$aspect))), 
                           spec_s1_co=rep(mean(brtdata$spec_s1_co),length(unique(brtdata$aspect))),
                           slope_per=rep(mean(brtdata$slope_per),length(unique(brtdata$aspect))),
                           spec_s2_co=rep(mean(brtdata$spec_s2_co),length(unique(brtdata$aspect))), 
                           DetOffset=aspoffset)
## elevation
elevvec <- seq(round(min(brtdata$elev), digits=0), round(max(brtdata$elev), digits=0), by=((round(max(brtdata$elev), digits=0) - round(min(brtdata$elev), digits=0))/20))
elevstep <- elevvec[2] - elevvec[1]
elevoffset <- vector()
for (i in 1:length(elevvec)){
  elevoffset[i] <- mean(brtdata$DetOffset[which((brtdata$elev>(elevvec[i]-elevstep/2)) & (brtdata$elev<(elevvec[i]+elevstep/2)))])
  if (is.nan(elevoffset[i])){    elevoffset[i] <- mean(brtdata$DetOffset)   }
  }
newdata.elev <- data.frame(pa_comm_code=rep("CS",length(elevvec)), 
                             aspect=rep("NE",length(elevvec)), 
                             elev=elevvec, 
                             spec_s1_co=rep(mean(brtdata$spec_s1_co),length(elevvec)),
                             slope_per=rep(mean(brtdata$slope_per),length(elevvec)),
                             spec_s2_co=rep(mean(brtdata$spec_s2_co),length(elevvec)), 
                             DetOffset=elevoffset)
## spec_s1_co
spec_s1_covec <- seq(round(min(brtdata$spec_s1_co), digits=0), round(max(brtdata$spec_s1_co), digits=0), by=((round(max(brtdata$spec_s1_co), digits=0) - round(min(brtdata$spec_s1_co), digits=0))/20))
spec_s1_costep <- spec_s1_covec[2] - spec_s1_covec[1]
spec_s1_cooffset <- vector()
for (i in 1:length(spec_s1_covec)){
  spec_s1_cooffset[i] <- mean(brtdata$DetOffset[which((brtdata$spec_s1_co>(spec_s1_covec[i]-spec_s1_costep/2)) & (brtdata$spec_s1_co<(spec_s1_covec[i]+spec_s1_costep/2)))])
  if (is.nan(spec_s1_cooffset[i])){ spec_s1_cooffset[i] <- mean(brtdata$DetOffset) }
}
newdata.spec_s1_co <- data.frame(pa_comm_code=rep("CS",length(spec_s1_covec)), 
                           aspect=rep("NE",length(spec_s1_covec)), 
                           elev=rep(mean(brtdata$elev),length(spec_s1_covec)), 
                           spec_s1_co=spec_s1_covec,
                           slope_per=rep(mean(brtdata$slope_per),length(spec_s1_covec)),
                           spec_s2_co=rep(mean(brtdata$spec_s2_co),length(spec_s1_covec)), 
                           DetOffset=spec_s1_cooffset)
## slope_per
slope_pervec <- seq(round(min(brtdata$slope_per), digits=0), round(max(brtdata$slope_per), digits=0), by=((round(max(brtdata$slope_per), digits=0) - round(min(brtdata$slope_per), digits=0))/20))
slope_perstep <- slope_pervec[2] - slope_pervec[1]
slope_peroffset <- vector()
for (i in 1:length(slope_pervec)){
  slope_peroffset[i] <- mean(brtdata$DetOffset[which((brtdata$slope_per>(slope_pervec[i]-slope_perstep/2)) & (brtdata$slope_per<(slope_pervec[i]+slope_perstep/2)))])
  if (is.nan(slope_peroffset[i])){ slope_peroffset[i] <- mean(brtdata$DetOffset) }
}
newdata.slope_per <- data.frame(pa_comm_code=rep("CS",length(slope_pervec)), 
                               aspect=rep("NE",length(slope_pervec)), 
                               elev=rep(mean(brtdata$elev),length(slope_pervec)), 
                               spec_s1_co=rep(mean(brtdata$spec_s1_co),length(slope_pervec)),
                               slope_per=slope_pervec,
                               spec_s2_co=rep(mean(brtdata$spec_s2_co),length(slope_pervec)), 
                               DetOffset=slope_peroffset)
## spec_s2_co
spec_s2_covec <- seq(round(min(brtdata$spec_s2_co), digits=0), round(max(brtdata$spec_s2_co), digits=0), by=((round(max(brtdata$spec_s2_co), digits=0) - round(min(brtdata$spec_s2_co), digits=0))/20))
spec_s2_costep <- spec_s2_covec[2] - spec_s2_covec[1]
spec_s2_cooffset <- vector()
for (i in 1:length(spec_s2_covec)){
  spec_s2_cooffset[i] <- mean(brtdata$DetOffset[which((brtdata$spec_s2_co>(spec_s2_covec[i]-spec_s2_costep/2)) & (brtdata$spec_s2_co<(spec_s2_covec[i]+spec_s2_costep/2)))])
  if (is.nan(spec_s2_cooffset[i])){ spec_s2_cooffset[i] <- mean(brtdata$DetOffset) }
}
newdata.spec_s2_co <- data.frame(pa_comm_code=rep("CS",length(spec_s2_covec)), 
                                aspect=rep("NE",length(spec_s2_covec)), 
                                elev=rep(mean(brtdata$elev),length(spec_s2_covec)), 
                                spec_s1_co=rep(mean(brtdata$spec_s1_co),length(spec_s2_covec)),
                                slope_per=rep(mean(brtdata$slope_per),length(spec_s2_covec)),
                                spec_s2_co=spec_s2_covec, 
                                DetOffset=spec_s2_cooffset)

# bootstrap brt models to get SE around each predictor in the simplified model for plotting
bootpreds.pacc <- matrix(nrow=100, ncol=nrow(newdata.pacc), byrow=TRUE)
bootpreds.aspect <- matrix(nrow=100, ncol=nrow(newdata.aspect), byrow=TRUE)
bootpreds.elev <- matrix(nrow=100, ncol=nrow(newdata.elev), byrow=TRUE)
bootpreds.spec_s1_co <- matrix(nrow=100, ncol=nrow(newdata.spec_s1_co), byrow=TRUE)
bootpreds.slope_per <- matrix(nrow=100, ncol=nrow(newdata.slope_per), byrow=TRUE)
bootpreds.spec_s2_co <- matrix(nrow=100, ncol=nrow(newdata.spec_s2_co), byrow=TRUE)
# repeat for each predictor in simplified model

# this produces 100 bootstraps. NOTE: THIS WILL TAKE A LONG TIME TO RUN
i <- 1
if (i<101){
  for (i in 1:100){
    bootdat <- brtdata[sample(1:nrow(brtdata), replace = TRUE),]
    brt.best.b <- try(gbm.step(data=bootdat, gbm.x = spec.varno, gbm.y = which(colnames(brtdata)=="TotSing"), offset=brtdata$DetOffset,
                               family="poisson", fold.vector=brtdata$LongFold, tree.complexity = spec.simp$interaction.depth,
                               learning.rate = 0.01, bag.fraction = spec.simp$bag.fraction), silent=T)
    if (class(brt.best.b)=="gbm"){
      bootpreds.pacc[i,] <- predict.gbm(brt.best.b, newdata.pacc[,1:6], n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.pacc$DetOffset)
      bootpreds.aspect[i,] <- predict.gbm(brt.best.b, newdata.aspect, n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.aspect$DetOffset)  
      bootpreds.elev[i,] <- predict.gbm(brt.best.b, newdata.elev, n.trees=brt.best.b$gbm.call$best.trees, type="response") *exp(newdata.elev$DetOffset) 
      bootpreds.spec_s1_co[i,] <- predict.gbm(brt.best.b, newdata.spec_s1_co, n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.spec_s1_co$DetOffset)  
      bootpreds.slope_per[i,] <- predict.gbm(brt.best.b, newdata.slope_per, n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.slope_per$DetOffset) 
      bootpreds.spec_s2_co[i,] <- predict.gbm(brt.best.b, newdata.spec_s2_co, n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.spec_s2_co$DetOffset)  
      # repeat for each predictor in simplified model
    } else { # sometimes the models don't converge. This tells R to rerun this iteration with a new bootstrapped dataset
      i <- i-1
    }
    
  }
}

getse.boot <- function(x){
  sd(as.vector(x), na.rm=TRUE) # se of the original sample = sd of the bootstrapped sample
}


# produce marginal plots of Shrub from brt models
bootpreds.pacc.se <- sapply(as.data.frame(bootpreds.pacc), getse.boot)
bootpreds.aspect.se <- sapply(as.data.frame(bootpreds.aspect), getse.boot)
bootpreds.elev.se <- sapply(as.data.frame(bootpreds.elev), getse.boot)  
bootpreds.spec_s1_co.se <- sapply(as.data.frame(bootpreds.spec_s1_co), getse.boot) 
bootpreds.slope_per.se <- sapply(as.data.frame(bootpreds.slope_per), getse.boot) 
bootpreds.spec_s2_co.se <- sapply(as.data.frame(bootpreds.spec_s2_co), getse.boot)
# repeat for each predictor in simplified model

bootpreds.allbrt.se <- rbind(cbind("se.pacc",bootpreds.pacc.se), cbind("se.aspect",bootpreds.aspect.se),
                             cbind("se.elev",bootpreds.elev.se), cbind("se.spec_s1_co",bootpreds.spec_s1_co.se), 
                             cbind("se.slope_per",bootpreds.slope_per.se), cbind("se.spec_s2_co",bootpreds.spec_s2_co.se))
write.csv(bootpreds.allbrt.se, file=paste(pathtofiles, "BTBW_Bootstrapped_SE_brt.csv", sep=""))



############################################################################
## PRODUCE RESPONSE PLOTS FOR EACH PREDICTOR IN MODEL, WITH BOOTSTRAPPED 95% CONFIDENCE INTERVALS

# Produce the plot. This code is for continuous predictors (e.g., elevation, cover)
jpeg(file=paste(pathtofiles,"BTBW_Elevation.jpeg", sep=""), width =5, height = 4, units = 'in', res = 300)
par(pin=c(5,5), cex=1, cex.lab=1.2, mar=c(4, 4, 3, 2)) 
predbrt.elev <- predict.gbm(brt.best.b, newdata.elev, n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.elev$DetOffset)
predbrt.elev.L95 <- predbrt.elev - bootpreds.elev.se*1.96
predbrt.elev.L95[which(predbrt.elev.L95<0)] <- 0
plot(predbrt.elev ~ elevvec, xlab="Elevation", ylab="BTBW density", main="", bty="o", ylim=c(max(-0.2,round(min((predbrt.elev - bootpreds.elev.se*1.96)))),round(max((predbrt.elev + bootpreds.elev.se*1.96)))), type="l", lwd=2, col="black")
lines((predbrt.elev + bootpreds.elev.se*1.96)~elevvec, lty="dotted", lwd=2, col="grey50")
lines(predbrt.elev.L95~elevvec, lty="dotted", lwd=2, col="grey50")
dev.off()

# Produce the plot. This code is for categorical predictors (e.g., PA community type, aspect)
jpeg(file=paste(pathtofiles,"BTBW_PACommType.jpeg", sep=""), width =5, height = 4, units = 'in', res = 300)
par(pin=c(5,5), cex=1, cex.lab=1.2, mar=c(4, 4, 3, 2)) 
predbrt.pacc <- predict.gbm(brt.best.b, newdata.pacc, n.trees=brt.best.b$gbm.call$best.trees, type="response")*exp(newdata.pacc$DetOffset)
predbrt.pacc.L95 <- bootpreds.pacc.se*1.96
for (i in 1:length(predbrt.pacc.L95)){
  if (predbrt.pacc.L95[i] < predbrt.pacc[i]){ predbrt.pacc.L95[i] <- predbrt.pacc[i] }
}
paccx <- sort(unique(brtdata$pa_comm_code))
paccx.int <- seq(1, length(paccx), by=1)
plot(predbrt.pacc ~ paccx.int, xaxt="n", xlab="PA Community Type", ylab="BTBW density", main="", bty="o", ylim=c(max(-0.2,round(min((predbrt.pacc - bootpreds.pacc.se*1.96)))),round(max((predbrt.pacc + bootpreds.pacc.se*1.96)))), type="p", pch=19, col="black")
arrows(paccx.int, (predbrt.pacc + bootpreds.pacc.se*1.96), paccx.int, (predbrt.pacc -predbrt.pacc.L95), angle=90, code=3, length=0, col="grey50")
axis(side=1, at=paccx.int, labels=paccx) 
dev.off()

