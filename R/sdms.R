# =================================================================================================================
#
# Iberian Conservation
# Species Distribution Models
#
# =================================================================================================================

# with some code snippets borrowed from https://damariszurell.github.io/HU-GCIB/

# -----------------------------------------------------------------------------------------------------------------
# IMPORT and DEFINE
# -----------------------------------------------------------------------------------------------------------------
library(raster)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ecospat)
library(mecofun)
library(dismo)
library(gam)
library(GRaF)

wd <- "C:/Users/Leon/Google Drive/03_LSNRS/Projects/Iberian_Conservation/iberian_conservation"
setwd(wd)

species <- "ursusarctos"
#species <- "lynxpardinus"

# -----------------------------------------------------------------------------------------------------------------
# FUNCTIONS
# -----------------------------------------------------------------------------------------------------------------
# AUTOMATED PREDICTIONS (altered for weighted predictions)
make.preds <- function(model, newdata) {
  require(dismo)
  require(gam)
  require(rpart)
  require(randomForest)
  require(gbm)
  require(maxnet)
  require(kernlab)
    
  switch(class(model)[1],
         Bioclim = predict(model, newdata),
         Domain = predict(model, newdata),
         glm = predict(model, newdata, type='response', weights = weight),
         Gam = predict(model, newdata, type='response', weights = weight),
         rpart = predict(model, newdata),
         randomForest = predict(model, newdata, type='prob', na.action = na.omit)[,2],
         gbm = predict.gbm(model, newdata, 
                           n.trees=model$gbm.call$best.trees, type="response"),
         maxnet = predict(model, newdata, type="logistic"),
         gausspr = predict(model, newdata, type="probabilities"),
         graf = as.numeric(predict(model, newdata, type="response")[,1]))
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
  require(kernlab)
  require(GRaF)
  
  # Make k-fold data partitions
  ks <- kfold(X_train, k = kfold, by = X_train[,y_name])
  
  cross_val_preds <- data.frame(row = row.names(X_train), 
                                cross_val_preds = numeric(length = nrow(X_train))) 
  
  print(paste("Model: ", class(model)[1]))
  
  for (i in seq_len(kfold)) {
    cv_train <- X_train[ks != i,]
    cv_test <- X_train[ks == i,]
    
    # force equal presence-absence for ML algorithms
    if (class(model)[1]=='randomForest' | class(model)[1]=='rpart' | class(model)[1]=='gbm' | class(model)[1]=='gausspr' | class(model)[1]=='graf') {
      n_pre_cv <- dim(cv_train[cv_train$species == 1,])[1]
      cv_train_abs <- sample_n(cv_train[cv_train$species == 0,], n_pre_cv)
      cv_train <- rbind(cv_train_abs, cv_train[cv_train$species == 1,])
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
                     maxnet = maxnet(p=cv_train[,y_name], data = cv_train[,x_name]),
                     #gausspr = update(model, data=cv_train))
                     gausspr = kernlab::gausspr(species~., data = cv_train[,c(x_name, y_name)], kernel = "rbfdot"),
                     graf = GRaF::graf(y = cv_train[,y_name], x = cv_train[,x_name], method = "Laplace", opt.l = T))
    
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


# -----------------------------------------------------------------------------------------------------------------
# PREDICTOR VARIABLES
# -----------------------------------------------------------------------------------------------------------------

# Mask
mask <- raster("Data/Mask/EUROPE_MASK_10km.tif")
if(species_name == "Lynx pardinus"){
  mask <- raster("Data/Mask/IBERIA_MASK_10km.tif")
}

# CHELSA Bioclim
l_bioclim <- list.files("Data/Predictors/CHELSA_Bioclim", pattern = ".tif$", full.names = T)
bioclim <- raster::stack(l_bioclim)
names(bioclim) <- paste0("bio", 1:dim(bioclim)[3])

# LC fractions
l_lcf <- list.files("Data/Predictors/LC_Fractions", pattern = ".tif$", full.names = T)
lcf <- raster::stack(l_lcf)
names(lcf) <- c("f_artifical", "f_bare", "f_cropland", "f_highveg", "f_lowveg")

# DEM
dem <- raster("Data/Predictors/DEM/GTOPO30_10km.tif_EUROPE.tif")
names(dem) <- "dem"

# Landsat Spectral Temporal Metrics
l_stm <- list.files("Data/Predictors/Landsat_STMs", pattern = ".tif$", full.names = T)
stm <- raster::stack(l_stm)
names(stm) <- apply(expand.grid(c("tcb", "tcg", "tcw"), c("sos", "pos", "eos")), 1, paste, collapse=".")


brick_preds <- brick(list(bioclim, lcf, dem, stm))

if(species_name == "Lynx pardinus"){
  brick_preds <- crop(brick_preds, mask)
}

brick_preds <- raster::mask(brick_preds, mask, maskvalue = 0)

# standardise data (z-transformation)
brick_preds_val <- getValues(brick_preds)
brick_preds_val <- scale(brick_preds_val)
brick_preds <- setValues(brick_preds, brick_preds_val)


# -----------------------------------------------------------------------------------------------------------------
# SPECIES DATA
# -----------------------------------------------------------------------------------------------------------------

# read presence/absence rasters
r_pre <- raster(paste0("Data/GBIF/", species, "_europe_1990-2020_presence_thinned.tif"))
r_abs <- raster(paste0("Data/GBIF/", species, "_europe_1990-2020_absence_thinned.tif"))

# mask absences closer than specified buffer to presences
buffer <- raster(paste0("Data/GBIF/", species, "_europe_1990-2020_presence_thinned_buffer.tif"))
r_abs <- mask(r_abs, buffer, maskvalue=1)
r_pre[r_pre == 1] <- 2

# convert to presence (1) - absence (0) - non (NA)
pre_abs <- brick(c(r_pre, r_abs))
pre_abs <- as.matrix(pre_abs)
pre_abs <- apply(pre_abs, 1, FUN=sum, na.rm = T)
pre_abs[pre_abs == 0] <- NA
pre_abs <- pre_abs-1

# cast predictor raster to dataframe and add species presence (1/0)
df <- as.data.frame(brick_preds, xy = TRUE)
df[, "species"] <- pre_abs   
df <- drop_na(df)

# specify predictor and response variable
X_names <- names(df)[4:length(df)-1]
y_name <- "species"

# variable correlations
#var_names <- append(X_names, y_name)
#cor_mat <- cor(df[var_names], method='spearman')
#corrplot(cor_mat, tl.srt=45, tl.col="black")

# apply variable selection
var_sel <- mecofun::select07_cv(X = df[, X_names], y = df[, y_name], kfold = 5, threshold = 0.7)
pred_sel <- var_sel$pred_sel
print(pred_sel)

# plot correlations among selected predictors
cor_mat <- cor(df[append(pred_sel, y_name)], method='spearman')
corrplot(cor_mat, tl.srt=45, tl.col="black")

# check how many data points are needed given the number of predictors
print(dim(df[df$species == 1,]))  # returns 867
print(length(pred_sel) * 10 * 2)  # returns 260


# -----------------------------------------------------------------------------------------------------------------
# MODEL FITTING
# -----------------------------------------------------------------------------------------------------------------

# specify number of absence points 
n_pre <- dim(df[df$species == 1,])[1]

if (n_pre * 10 < dim(df[df$species == 0,])[1]) {
  n_abs <- n_pre * 10
} else {
  n_abs <- dim(df[df$species == 0,])[1]
}

# weighting presence-absence data
abs_weights <- n_pre / n_abs
df_sub <- sample_n(df[df$species == 0,], n_abs)
df_sub <- rbind(df_sub, df[df$species == 1,])

df_sub$weight <- ifelse(df_sub$species == 1, 1, abs_weights)  # abs_weight or 1 if equally 


# (1.1) Generalized Linear Models (GLMs)
m_glm <- step(
  glm(as.formula(paste('species ~', paste(pred_sel, paste0('+ I(', pred_sel,'^2)'), collapse=' + '))),
      family='binomial', data = df_sub, weights = weight)
)
saveRDS(m_glm, paste0("R/models/", species, "_m_glm.rds"))


# (1.2) Generalized Additive Models (GAMs)
m_gam <- gam( as.formula(
  paste('species ~', paste(paste0('s(',pred_sel,',df=4)'), collapse=' + '))),
  family='binomial', data=df_sub, weights = weight)
saveRDS(m_gam, paste0("R/models/", species, "_m_gam.rds"))


# (1.3) BIOCLIM
m_bc <- bioclim(brick_preds[[pred_sel]], df_sub[df_sub$species == 1, c('x','y')])
saveRDS(m_bc, paste0("R/models/", species, "_m_bc.rds"))


# selected machine learners require equal sample size of presence-absence
#df_sub_ml <- sample_n(df[df$species == 0,], n_pre)
#df_sub_ml <- rbind(df_sub_ml, df[df$species == 1,])
 
# according to Barbet-Massin et al. (2012): "Same as number of presences, 10 runs when less than 1000 PA with an equal weight for presences and absences"
c_brt <- c()
c_gp <- c()

for (i in 1:10) {
  df_sub_i <- sample_n(df_sub[df_sub$species == 0,], n_pre)
  df_sub_i <- rbind(df_sub_i, df_sub[df_sub$species == 1,])
  
  # (2.1) Boosted regression trees (BRT)
  m_brt <- gbm.step(data = df_sub_i, 
                    gbm.x = pred_sel,
                    gbm.y = 'species', 
                    family = 'bernoulli',
                    tree.complexity = 2,
                    bag.fraction = 0.75,
                    learning.rate = 0.001,
                    verbose = F)
  assign(paste0("m_brt", i), m_brt)
  c_brt[i] <- paste0("m_brt", i)
  
  # (2.2) Gaussian Process Classification (GPC)
  m_gp <- GRaF::graf(y = df_sub_i$species, x = df_sub_i[,pred_sel], opt.l = TRUE, method = "Laplace")
  assign(paste0("m_gp", i), m_gp)
  c_gp[i] <- paste0("m_gp", i)
  
  saveRDS(m_brt, paste0("R/models/", species, "_m_brt", i,".rds"))
  saveRDS(m_gp, paste0("R/models/", species, "_m_gp", i,".rds"))
}



# slow, but optimises lengthscale of RBF kernel
#m_gp <- GRaF::graf(y = df_sub_ml$species, x = df_sub_ml[,pred_sel], method = "Laplace")
# fast, but only approximates lenghtscale by input features

#gpc_pred_prob <- raster::predict(m_gpc1, brick_preds[[pred_sel]], type="response", progress = "text")  # probability scale predictions
#gpc_pred_lat <- raster::predict(m_gpc1, brick_preds[[pred_sel]], type="latent", CI = 'std', progress = "text")  # mean & std on Gaussian scale 
#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#plot(gpc_pred_prob[[1]], col=jet.colors(100))
#plot(gpc_pred_lat[[2]]^2, col=jet.colors(100))  # 'uncertainty', i.e. std of the multivariate Gaussian around the mean prediction

# ----- OLD MODELS
# (4) Domain
#m_dom <- domain(brick_preds[[pred_sel]], df_sub[df_sub$species == 1, c('x','y')])

# (5) CARTs
#library(rpart)
#m_cart <- rpart( as.formula(
#paste('species ~', paste(pred_sel, collapse=' + '))),
#data=df_sub_ml, control=rpart.control(minsplit=20,xval=10))


# (6) Random Forest (RF)
#library(randomForest)
#m_rf <- randomForest(x=df_sub_ml[,pred_sel], y=as.factor(df_sub_ml$species), 
#ntree=1000, importance =T, na.action = na.omit)
#graf = predict(m_gp, df_sub[, pred_sel], type="response")



# -----------------------------------------------------------------------------------------------------------------
# MODEL ASSESSMENT
# -----------------------------------------------------------------------------------------------------------------
all_models <- c('m_glm', 'm_gam', 'm_bc', 'm_brt', 'm_gp')
stat_models <- c('m_glm', 'm_gam', 'm_bc')
ml_models <- c(c_brt, c_gp)


# (1) internal (training data) performance
# make predictions
train_preds_stat <- sapply(stat_models, FUN=function(m){make.preds(eval(parse(text=m)), df_sub[, pred_sel])})
train_preds_ml <- sapply(ml_models, FUN=function(m){make.preds(eval(parse(text=m)), df_sub[, pred_sel])})

# average ml predictions
brt_means <- train_preds_ml %>%
  as.data.frame() %>%
  select(c(1,1:10)) %>%
  rowMeans()
gp_means <- train_preds_ml %>%
  as.data.frame() %>%
  select(c(1,11:20)) %>%
  rowMeans()

# combine all models
train_preds <- data.frame(cbind(train_preds_stat, brt_means, gp_means))
names(train_preds) <- all_models

# assess performance
train_perf <- data.frame(sapply(all_models, FUN=function(x){calc.eval(df_sub, 'species', train_preds[,x])}))
print(train_perf)
write.csv2(apply(train_perf, 2, as.character), paste("Data/SDMs/", species, "_internal_model_validation_sel07-th05.csv"))


# (2) cross-validated performance
# make predictions
crossval_preds <- sapply(all_models, FUN = function(x) {crossval.preds(model = eval(parse(text=x)),
                                                                     X_train = df_sub, 
                                                                     y_name = 'species', 
                                                                     x_name = pred_sel, 
                                                                     X_raster = brick_preds, 
                                                                     colname_coord = c('x','y'), 
                                                                     kfold=5)})
# assess performance
crossval_perf <- data.frame(sapply(all_models, FUN=function(x){calc.eval(df_sub, 'species', crossval_preds[,x])}))
print(crossval_perf)
write.csv2(apply(crossval_perf, 2, as.character), paste("Data/SDMs/", species, "_k5-fold_crossval_model_validation_sel07-th05.csv"))


# -----------------------------------------------------------------------------------------------------------------
# MODEL ENSEMBLES
# -----------------------------------------------------------------------------------------------------------------
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
ensemble_perf <- sapply(names(ensemble_preds)[1:4], FUN=function(x){calc.eval(df_sub, 'species', ensemble_preds[,x])})
summary(ensemble_preds)
print(ensemble_perf)
write.csv2(apply(ensemble_perf,2,as.character), paste("Data/SDMs/", species, "_ensemble_model_validation_sel07-th05.csv"))


# ----------------------------------------------------------------------------------------------------
# MAP PREDICTIONS
# ----------------------------------------------------------------------------------------------------
EU <- data.frame(rasterToPoints(brick_preds[[pred_sel]]))
EU[is.na(EU)] <- 0
EU <- drop_na(EU)

# make predictions
img_pred <- data.frame(EU[,1:2], sapply(stat_models, FUN=function(m){make.preds(eval(parse(text=m)), EU[, pred_sel])}))
img_pred_ml <- data.frame(EU[,1:2], sapply(ml_models, FUN=function(m){make.preds(eval(parse(text=m)), EU[, pred_sel])}))
img_pred_brt <- img_pred_ml %>%
  select(c_brt) %>%
  rowMeans()
img_pred_gp <- img_pred_ml %>%
  select(c_gp) %>%
  rowMeans()
img_pred <- cbind(img_pred, img_pred_brt, img_pred_gp)
names(img_pred) <- c("x", "y", all_models)

# Binarise predictions of all algorithms
env_preds_bin <- data.frame(EU[,1:2], sapply(all_models, FUN=function(x){ifelse(img_pred[,x]>= unlist(crossval_perf['thresh',x]),1,0)}))

# Make rasters from predictions:
r_preds <- rasterFromXYZ(img_pred, crs=crs(brick_preds))
r_preds_bin <- rasterFromXYZ(env_preds_bin, crs=crs(brick_preds))

writeRaster(r_preds, paste0("Data/SDMs/", species, "_SDM_probability_models.tif"), "GTiff", overwrite = T)
writeRaster(r_preds_bin, paste0("Data/SDMs/", species, "_SDM_binary_models.tif"), "GTiff", overwrite = T)

# We make ensembles:    
env_ensemble <- data.frame(EU[,1:2], make.ensemble(img_pred[,-c(1:2)], unlist(crossval_perf['TSS',]), unlist(crossval_perf['thresh',])))

# Make rasters from ensemble predictions:
r_ens <- rasterFromXYZ(env_ensemble, crs=crs(brick_preds))
writeRaster(r_ens, paste0("Data/SDMs/", species, "_SDM_ensemble_mn_md_wmn_cav_sdp.tif"), "GTiff", overwrite = T)

# Binarise ensemble predictions
env_ensemble_bin <- data.frame(EU[,1:2], sapply(c('mean_prob', 'median_prob', 'wmean_prob'), FUN=function(x){ifelse(env_ensemble[,x]>= unlist(ensemble_perf['thresh',x]),1,0)}))

# Make rasters:
r_ens_bin <- rasterFromXYZ(env_ensemble_bin, crs=crs(brick_preds))
writeRaster(r_ens_bin, paste0("Data/SDMs/", species, "_SDM_ensembles_binary_mn_md_wmn.tif"), "GTiff", overwrite = T)

# quick visual checks (just temporarily to keep)
prob_models <- stack("Data/SDMs/lynxpardinus_SDM_probability_models.tif")
bin_models <- stack("Data/SDMs/lynxpardinus_SDM_binary_models.tif")
names(prob_models) <- c("GLM", "GAM", "Bioclim", "BRT", "GP")
names(bin_models) <- c("GLM", "GAM", "Bioclim", "BRT", "GP")
spplot(prob_models)
spplot(bin_models)


