# ====================================================================================================
#
# 2 GLOBAL CHANGE IMPACTS ON BIODIVERSITY:
# Species Distribution Models
#
# ====================================================================================================

# with some code snippets borrowed from https://damariszurell.github.io/HU-GCIB/

# ----------------------------------------------------------------------------------------------------
# IMPORT and DEFINE
# ----------------------------------------------------------------------------------------------------
library(raster)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ecospat)


wd <- "/Users/leonnill/Google Drive/01_MSc_GCG/MSc6_GCIB/Project/Data"
setwd(wd)


# ----------------------------------------------------------------------------------------------------
# FUNCTIONS
# ----------------------------------------------------------------------------------------------------

# VARIABLE SELECTION AND MULTICOLLINEARITY
select07 <- function(pred_names, response_name, data, cor_mat=NULL, threshold=0.7){
  
  # Function for calculating AIC - we use univariate GLMs with linear and quadratic terms
  var.imp <- function (predictor, response)
  {
    AIC(glm(response ~ predictor + I(predictor^2), binomial))
  }
  
  # Calculate AIC for all predictor variables
  aic_imp <- apply(data[pred_names], 2, var.imp, response= data[,response_name])
  
  # Names of sorted variables
  sort_imp <- names(sort(aic_imp))
  
  # Calculate correlation matrix if not provided in function call
  if (is.null(cor_mat)) {
    cor_mat <- cor(data[pred_names], method='spearman')
  }
  
  # Identifies correlated variable pairs:
  diag(cor_mat)=NA
  pairs <- which(abs(cor_mat)>= threshold, arr.ind=T) 
  
  # Identify which variables should be excluded
  exclude <- NULL
  for (i in 1:length(sort_imp))
  {
    if ((sort_imp[i] %in% row.names(pairs))&
        ((sort_imp[i] %in% exclude)==F)) {
      cv <- cor_mat[setdiff(row.names(cor_mat),exclude),sort_imp[i]]
      cv <- cv[setdiff(names(cv),sort_imp[1:i])]
      exclude <- c(exclude,names(which((abs(cv)>=threshold)))) 
    }
  }
  
  # Select set of weakly correlated predictors:
  pred_sel <- sort_imp[!(sort_imp %in% exclude)]
  
  # Return list with AIC, correlation matrix, and final predictors:
  return(list(AIC=sort(aic_imp), cor_mat=cor_mat, pred_sel=pred_sel))
}


# AUTOMATED PREDICTIONS (altered for weighted predictions)
make.preds <- function(model, newdata) {
  require(dismo)
  require(gam)
  require(rpart)
  require(randomForest)
  require(gbm)
  require(maxnet)
  
  switch(class(model)[1],
         Bioclim = predict(model, newdata),
         Domain = predict(model, newdata),
         glm = predict(model, newdata, type='response', weights = weight),
         Gam = predict(model, newdata, type='response', weights = weight),
         rpart = predict(model, newdata),
         randomForest = predict(model, newdata, type='prob', na.action = na.omit)[,2],
         gbm = predict.gbm(model, newdata, 
                           n.trees=model$gbm.call$best.trees, type="response"),
         maxnet = predict(model, newdata, type="logistic"))
}


# K-FOLD CROSSVALIDATION
crossval.preds <- function(model, X_train, y_name, x_name,
                           X_raster, colname_coord, kfold) {
  require(dismo)  
  require(gam)
  require(rpart)
  require(randomForest)
  require(gbm)
  require(maxnet)
  
  # Make k-fold data partitions
  ks <- kfold(X_train, k = kfold, by = X_train[,y_name])
  
  cross_val_preds <- data.frame(row = row.names(X_train), 
                                cross_val_preds = numeric(length = nrow(X_train))) 
  
  for (i in seq_len(kfold)) {
    cv_train <- X_train[ks != i,]
    cv_test <- X_train[ks == i,]
    
    # force equal presence-absence for ML algorithms
    if (class(model)[1]=='randomForest' | class(model)[1]=='rpart' | class(model)[1]=='gbm') {
      n_pre_cv <- dim(cv_train[cv_train$ursus_arctos == 1,])[1]
      cv_train_abs <- sample_n(cv_train[cv_train$ursus_arctos == 0,], n_pre_cv)
      cv_train <- rbind(cv_train_abs, cv_train[cv_train$ursus_arctos == 1,])
    }
    
    # Because we used the gbm.step() for BRTs, we need a small work-around:
    if (class(model)[1]=='gbm') {
      cv_train_gbm <- cv_train;
      names(cv_train_gbm)[names(cv_train_gbm)==y_name] <- 
        model$response.name
    }
    
    # We update the model for the new training data
    modtmp <- switch(class(model)[1],
                     Bioclim = bioclim(X_raster[[x_name]], cv_train[cv_train[, y_name]==1, colname_coord]),
                     Domain = domain(X_raster[[x_name]], cv_train[cv_train[, y_name]==1, colname_coord]),
                     glm = update(model, data=cv_train, weights = weight),
                     Gam = update(model, data=cv_train, weights = weight),
                     rpart = update(model, data=cv_train),
                     randomForest = update(model, data=cv_train),                
                     gbm = gbm(model$call, 'bernoulli', data=cv_train_gbm[,c(x_name,model$response.name)], n.trees=model$gbm.call$best.trees, shrinkage=model$gbm.call$learning.rate, bag.fraction=model$gbm.call$bag.fraction, interaction.depth=model$gbm.call$tree.complexity),
                     maxnet = maxnet(p= cv_train[,y_name], data= cv_train[,x_name]))
    
    # We make predictions for k-fold:
    if (class(model)[1]=='gbm') {
      cross_val_preds[which(ks==i),2] <- 
        predict.gbm(modtmp, cv_test[, x_name], n.trees=model$gbm.call$best.trees, type="response")
    } else {
      cross_val_preds[which(ks==i),2] <- make.preds(modtmp, cv_test[, x_name])
    }
  }
  return(cross_val_preds[,2])
  #return(cross_val_preds)
}


# MODEL PERFORMANCE MEASURES
calc.eval <- function(dat, colname_species, preds, thresh_method='MaxSens+Spec'){
  require(PresenceAbsence)
  require(dismo)  
  
  # Helper functions.
  ## True Skill Statistic:
  TSS = function(cmx){
    PresenceAbsence::sensitivity(cmx, st.dev=F) + 
      PresenceAbsence::specificity(cmx, st.dev=F) - 1
  }
  
  ## Explained deviance - the function calc.deviance() is taken from the dismo package:
  d.square <- function(obs, pred, family='binomial'){
    pred <- ifelse(pred<.00001,.00001,ifelse(pred>.9999,.9999,pred))
    
    null_pred <- rep(mean(obs), length(obs))
    
    1 - (calc.deviance(obs, pred, family=family) / 
           calc.deviance(obs, null_pred, family=family))
  }
  
  # Prepare data set to optimise threshold for binarising:
  thresh_dat <- data.frame(ID=seq_len(nrow(dat)), 
                           obs = dat[, colname_species],
                           pred = preds)
  
  # Find optimal threshold using the package PresenceAbsence
  thresh <- optimal.thresholds(DATA= thresh_dat)
  
  # Prepare confusion matrix
  cmx_maxSSS <- cmx(DATA= thresh_dat, threshold=thresh[thresh$Method==thresh_method,2])
  
  # Generate output data frame with performance statistics and optimal threshold:
  data.frame(AUC = PresenceAbsence::auc(thresh_dat, st.dev=F),
             TSS = TSS(cmx_maxSSS), 
             Sens = PresenceAbsence::sensitivity(cmx_maxSSS, st.dev=F),
             Spec = PresenceAbsence::specificity(cmx_maxSSS, st.dev=F),
             PCC = PresenceAbsence::pcc(cmx_maxSSS, st.dev=F), 
             D2 = d.square(thresh_dat$obs, thresh_dat$pred),
             thresh = thresh[thresh$Method==thresh_method,2])
}


# ----------------------------------------------------------------------------------------------------
# PREDICTOR VARIABLES
# ----------------------------------------------------------------------------------------------------

# Mask
mask <- raster("Mask/EUROPE_MASK_10km.tif")

# CHELSA Bioclim
l_bioclim <- list.files("Predictor Variables/CHELSA_Bioclim", pattern = ".tif$", full.names = T)
bioclim <- raster::stack(l_bioclim)
names(bioclim) <- paste0("bio", 1:dim(bioclim)[3])

# LC fractions
l_lcf <- list.files("Predictor Variables/LC_Fractions", pattern = ".tif$", full.names = T)
lcf <- raster::stack(l_lcf)
names(lcf) <- c("f_artifical", "f_bare", "f_cropland", "f_highveg", "f_lowveg")

# DEM
dem <- raster("Predictor Variables/DEM/GTOPO30_10km.tif_EUROPE.tif")
names(dem) <- "dem"

# Landsat Spectral Temporal Metrics
l_stm <- list.files("Predictor Variables/Landsat_STMs", pattern = ".tif$", full.names = T)
stm <- raster::stack(l_stm)
names(stm) <- apply(expand.grid(c("tcb", "tcg", "tcw"), c("sos", "pos", "eos")), 1, paste, collapse=".")


brick_preds <- brick(list(bioclim, lcf, dem, stm))
brick_preds <- raster::mask(brick_preds, mask, maskvalue = 0)

# standardise data (z-transformation)
brick_preds_val <- getValues(brick_preds)
brick_preds_val <- scale(brick_preds_val)
brick_preds <- setValues(brick_preds, brick_preds_val)


# ----------------------------------------------------------------------------------------------------
# SPECIES DATA
# ----------------------------------------------------------------------------------------------------

r_pre <- raster("GBIF/UrsusArctos_EUROPE_10km_presence_thinned.tif")
r_abs <- raster("GBIF/UrsusArctos_EUROPE_10km_absence_thinned.tif")

# with buffer
buffer <- raster("GBIF/UrsusArctos_EUROPE_10km_presence_thinned_BUFFER_40km.tif")
r_abs <- mask(r_abs, buffer, maskvalue=1)
r_pre[r_pre == 1] <- 2

pre_abs <- brick(c(r_pre, r_abs))
pre_abs <- as.matrix(pre_abs)
pre_abs <- apply(pre_abs, 1, FUN=sum, na.rm = T)
pre_abs[pre_abs == 0] <- NA
pre_abs <- pre_abs-1

df_preds <- as.data.frame(brick_preds, xy = TRUE)
df <- cbind(df_preds, "ursus_arctos" = pre_abs)
df <- drop_na(df)

X_names <- names(df)[4:length(df)-1]
y_name <- "ursus_arctos"

# variable correlations
var_names <- append(X_names, y_name)
cor_mat <- cor(df[var_names], method='spearman')
corrplot(cor_mat, tl.srt=45, tl.col="black")

var_sel <- select07(pred_names=X_names, response_name=y_name, data=df, cor_mat=cor_mat, threshold=0.5)
pred_sel <- var_sel$pred_sel
print(pred_sel)
cor_mat <- cor(df[append(pred_sel, y_name)], method='spearman')
corrplot(cor_mat, tl.srt=45, tl.col="black")

# Check how many predictors could potentially be included in the models
print(dim(df[df$ursus_arctos == 1,]))  # returns 481
print(length(pred_sel) * 10 * 2)  # returns 260


# ----------------------------------------------------------------------------------------------------
# MODEL FITTING
# ----------------------------------------------------------------------------------------------------

# weighting presence-absence data
n_pre <- dim(df[df$ursus_arctos == 1,])[1]
n_abs <- 5000
abs_weights <- n_pre / n_abs

df_sub <- sample_n(df[df$ursus_arctos == 0,], n_abs)
df_sub <- rbind(df_sub, df[df$ursus_arctos == 1,])

df_sub$weight <- ifelse(df_sub$ursus_arctos == 1, 1, abs_weights)  # abs_weight or 1 if equally 

# (1) Generalized Linear Models (GLMs)
m_glm <- step(
    glm(as.formula(paste('ursus_arctos ~', paste(pred_sel, paste0('+ I(', pred_sel,'^2)'), collapse=' + '))),
        family='binomial', data = df_sub, weights = weight)
  )


# (2) Generalized Additive Models (GAMs)
library(gam)
m_gam <- gam( as.formula(
  paste('ursus_arctos ~', paste(paste0('s(',pred_sel,',df=4)'), collapse=' + '))),
  family='binomial', data=df_sub, weights = weight)


# (3) BIOCLIM
library(dismo)
m_bc <- bioclim(brick_preds[[pred_sel]], df_sub[df_sub$ursus_arctos == 1, c('x','y')])


# (4) Domain
m_dom <- domain(brick_preds[[pred_sel]], df_sub[df_sub$ursus_arctos == 1, c('x','y')])


# Machine Learning Methods require equal sample size of presence-absence
n_pre <- dim(df[df$ursus_arctos == 1,])[1]
n_abs <- dim(df[df$ursus_arctos == 1,])[1]

df_sub_ml <- sample_n(df[df$ursus_arctos == 0,], n_abs)
df_sub_ml <- rbind(df_sub_ml, df[df$ursus_arctos == 1,])

# (5) CARTs
library(rpart)
m_cart <- rpart( as.formula(
  paste('ursus_arctos ~', paste(pred_sel, collapse=' + '))),
  data=df_sub_ml, control=rpart.control(minsplit=20,xval=10))


# (6) Random Forest (RF)
library(randomForest)
m_rf <- randomForest(x=df_sub_ml[,pred_sel], y=as.factor(df_sub_ml$ursus_arctos), 
                     ntree=1000, importance =T, na.action = na.omit)


# (7) Boosted regression trees (BRT)
library(gbm)
m_brt <- gbm.step(data = df_sub_ml, 
                  gbm.x = pred_sel,
                  gbm.y = 'ursus_arctos', 
                  family = 'bernoulli',
                  tree.complexity = 2,
                  bag.fraction = 0.75,
                  learning.rate = 0.001,
                  verbose=F)


# ----------------------------------------------------------------------------------------------------
# MODEL ASSESSMENT
# ----------------------------------------------------------------------------------------------------
c_models <- c('m_glm', 'm_bc', 'm_dom', 'm_gam', 'm_cart', 'm_rf', 'm_brt')

# MAKE PREDICTIONS
train_preds <- sapply(c_models, FUN=function(m){make.preds(eval(parse(text=m)), df_sub[, pred_sel])})

# INTERNAL MODEL PERFORMANCE
train_perf <- data.frame(sapply(c_models, FUN=function(x){calc.eval(df_sub, 'ursus_arctos', train_preds[,x])}))
print(train_perf)
write.csv2(apply(train_perf,2,as.character), "internal_model_validation_sel07-th05.csv")

# CROSS-VALIDATED PREDICTIONS
crossval_preds <- sapply(c_models, FUN = function(x) {crossval.preds(model = eval(parse(text=x)),
                                                                  X_train = df_sub, 
                                                                  y_name = 'ursus_arctos', 
                                                                  x_name = pred_sel, 
                                                                  X_raster = brick_preds, 
                                                                  colname_coord = c('x','y'), 
                                                                  kfold=5)})


# K-FOLD MODEL PERFORMANCE
crossval_perf <- data.frame(sapply(c_models, FUN=function(x){calc.eval(df_sub,
                                                                       'ursus_arctos',
                                                                       crossval_preds[,x])}))
print(crossval_perf)
write.csv2(apply(crossval_perf,2,as.character), "k5-fold_crossval_model_validation_sel07-th05.csv")


# ----------------------------------------------------------------------------------------------------
# MODEL ENSEMBLES
# ----------------------------------------------------------------------------------------------------
make.ensemble <- function(preds, eval_metric, thresh){
  # "preds" is a data.frame containing predictions for different algorithms.
  # "eval_metric" is a vector with same length as number of columns in preds. It provides the evaluation metric used for weighting probabilities.
  # "thresh" is a vector with same length as number of columns in preds. It provides the algorithm-specific threshold for making binary predictions.
  
  data.frame(mean_prob = rowMeans(preds),
             median_prob = apply(preds,1,median),
             wmean_prob = apply(preds,1,weighted.mean, w=eval_metric),
             committee_av = rowSums(sapply(seq_len(ncol(preds)), FUN=function(x){
               ifelse(preds[,x]>=thresh[x],1,0)}))/ncol(preds),
             sd_prob = apply(preds,1,sd))
}

# Make ensemble predictions:
ensemble_preds <- make.ensemble(crossval_preds,
                                unlist(crossval_perf['TSS',]),
                                unlist(crossval_perf['thresh',]))

# Evaluate ensemble predictions:
ensemble_perf <- sapply(names(ensemble_preds)[1:4], FUN=function(x){calc.eval(df_sub, 'ursus_arctos', ensemble_preds[,x])})
summary(ensemble_preds)
print(ensemble_perf)
write.csv2(apply(ensemble_perf,2,as.character), "ensemble_model_validation_sel07-th05.csv")


# ----------------------------------------------------------------------------------------------------
# MAP PREDICTIONS
# ----------------------------------------------------------------------------------------------------
EU <- data.frame(rasterToPoints(brick_preds[[pred_sel]]))
EU <- drop_na(EU)

# We make predictions of all models:
env_preds <- data.frame(EU[,1:2], sapply(c_models, FUN=function(m){make.preds(eval(parse(text=m)), EU)}))

# Binarise predictions of all algorithms
env_preds_bin <- data.frame(EU[,1:2], sapply(c_models, FUN=function(x){ifelse(env_preds[,x]>= unlist(crossval_perf['thresh',x]),1,0)}))

# Make rasters from predictions:
r_preds <- rasterFromXYZ(env_preds, crs=crs(brick_preds))
r_preds_bin <- rasterFromXYZ(env_preds_bin, crs=crs(brick_preds))

writeRaster(r_preds, "SDMs/SDM_probability_single.tif", "GTiff")
writeRaster(r_preds_bin, "SDMs/SDM_binary_single.tif", "GTiff")

# Map predicted occurrence probabilities:
spplot(r_preds)


# We make ensembles:    
env_ensemble <- data.frame(EU[,1:2], make.ensemble(env_preds[,-c(1:2)], unlist(crossval_perf['TSS',]), unlist(crossval_perf['thresh',])))

# Make rasters from ensemble predictions:
r_ens <- rasterFromXYZ(env_ensemble, crs=crs(brick_preds))
writeRaster(r_ens, "SDMs/SDM_ensembles_mn_md_wmn_cav_sdp.tif", "GTiff")

# Map continuous ensemble predictions:
spplot(r_ens[[1:4]])


# Binarise ensemble predictions
env_ensemble_bin <- data.frame(EU[,1:2], sapply(c('mean_prob', 'median_prob', 'wmean_prob'), FUN=function(x){ifelse(env_ensemble[,x]>= unlist(ensemble_perf['thresh',x]),1,0)}))

# Make rasters:
r_ens_bin <- rasterFromXYZ(env_ensemble_bin, crs=crs(brick_preds))
writeRaster(r_ens_bin, "SDMs/SDM_ensembles_binary_mn_md_wmn.tif", "GTiff")

# Map predicted presence from ensembles:
spplot(r_ens_bin)  


