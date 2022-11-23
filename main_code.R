library(tidyverse)   # tidy data
library(dlookr)      # for exploratory data analysis and imputation
library(visdat)      # for visualizing NAs
library(plotly)      # for interactive visualization
library(missRanger)  # to generate and impute NAs
library(GGally)      # vis
library(hydroGOF)    # metrics calculation
library(lime)       
library(gbm)        # model building
library(caret)      # model building
library(rsample)    # split the data into train and test set

#Working with DGWL data, this makes this work different from the paper
##First, step is processing the DGWL data and append it to master file
###Import the data


DGWL <- read_csv("Dataset/dgwlData.csv")
summary(DGWL)

#Convert wide to long and merge with the master dataframe
DGWL_longer <- DGWL %>% 
  pivot_longer(cols = Barry:SWS, names_to = "Well.Name", values_to = "DGWL_m")
summary(DGWL_longer)

#merge the data frames
##rename the columns
colnames(DGWL_longer) <- c("Month", "Year", "wellName", "DGWL_m")

dgwlDataset <- merge(allData, DGWL_longer, by = c("Month", "Year", "wellName"))
summary(dgwlDataset)

write_csv(dgwlDataset, "dgwlData.csv")

#First part of the model building, is impute the missing DGWL data. This was accomplished, by building individual model for each well.

dgwlDataset <- dgwlDataset[,-16]    #removes GWLA

dgwlDataList <- split(dgwlDataset, dgwlDataset$wellName) 

#starts here for Github

windows()
plot_na_pareto(dgwlDataset, only_na = TRUE)

summary(dgwlDataset)

#visualize
dgwlDataset %>%
  select(Discharge, TRMM, SoilMoisture, GRACE_TWSA, DGWL_m) %>%
  GGally::ggpairs()

#Impute K with the minimum value, because the wells are located in bedrock
dgwlDataset$Conductivity[is.na(dgwlDataset$Conductivity)] <- 6.184
summary(dgwlDataset)

#Impute DGWL values using rpart
#reference: https://yuzar-blog.netlify.app/posts/2021-03-04-how-to-impute-missing-values-in-r/
imputedDgwlDatarpart <- dgwlDataset %>% 
  group_by(wellName) %>% 
  imputate_na(DGWL_m, method = "rpart")

imputedDgwlDatamice <- dgwlDataset %>% 
  group_by(wellName) %>% 
  imputate_na(DGWL_m, method = "mice") #mice gives the best result

summary(imputedDgwlDatamice)

plot(imputedDgwlDatamice)
plot(imputedDgwlDatarpart)
head(imputedDgwlDatamice)

dgwlDataset$dgwl <- imputedDgwlDatamice
dgwlDataset <- dgwlDataset[, -c(3, 16)]
summary(dgwlDataset)

#Model building
##No scaling 
#Split the data
set.seed(15)
split <- initial_split(dgwlDataset, 0.8)
train_dataset <- training(split)
test_dataset <- testing(split)

windows()
train_dataset %>%
  select(Discharge, TRMM, SoilMoisture, GRACE_TWSA, dgwl) %>%
  GGally::ggpairs()

gbm.fit <- gbm(formula = dgwl~.,
               distribution = "gaussian",
               data = train_dataset,
               n.trees = 20000,
               interaction.depth = 1,
               shrinkage = 0.001,
               cv.folds = 10,
               n.cores = NULL, 
               verbose = FALSE)  

sqrt(min(gbm.fit$train.error))

sqrt(min(gbm.fit$cv.error))     #RMSE

gbm.perf(gbm.fit, method = "cv")
# find index for n trees with minimum CV error
min_MSE <- which.min(gbm.fit$cv.error)

# 2nd train GBM model
set.seed(123)
gbm.fit2 <- gbm(
  formula = dgwl~.,
  distribution = "gaussian",
  data = train_dataset,
  n.trees = 10000,
  interaction.depth = 3,
  shrinkage = 0.1,
  cv.folds = 5,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
)

gbm.perf(gbm.fit2, method = "cv")
min_MSE <- which.min(gbm.fit2$cv.error)      # find index for n trees with minimum CV error
sqrt(min(gbm.fit2$cv.error))     #RMSE
sqrt(min(gbm.fit2$train.error))

#Just checking the test fits.
pred <- predict(gbm.fit2, n.trees = gbm.fit2$n.trees, train_dataset)
pred_test <- predict(gbm.fit2, n.trees = gbm.fit2$n.trees, test_dataset)

RMSE(pred_test, test_dataset$dgwl)

plot(pred_test, test_dataset$dgwl)
abline(coef=c(0,1), col="red", lwd=2)

#Tune model
# create hyperparameter grid
hyper_grid <- expand.grid(shrinkage = c(.01, .1, .3),
                          interaction.depth = c(1, 3, 5),
                          n.minobsinnode = c(5, 10, 15),
                          bag.fraction = c(.65, .8, 1),
                          optimal_trees = 0,          # a place to dump results
                          min_RMSE = 0                # a place to dump results
)                      

# total number of combinations
nrow(hyper_grid)

# randomize data (timeseries data need to be randomized)
random_index <- sample(1:nrow(train_dataset), nrow(train_dataset))
random_gwla_train <- train_dataset[random_index, ]

# grid search 
for(i in 1:nrow(hyper_grid)) {
  # reproducibility
  set.seed(123)
  # train model
  gbm.tune <- gbm(
    formula = dgwl ~ .,
    distribution = "gaussian",
    data = random_gwla_train,
    n.trees = 5000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .75, #only uses 75% of the training data, the remaining 25% used for evaluation of performance
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

hyper_grid %>% 
  dplyr::arrange(min_RMSE) %>%
  head(10)

# create hyperparameter grid
hyper_grid2 <- expand.grid(shrinkage = c(.005, .01, .1),
                           interaction.depth = c(5, 7, 10),
                           n.minobsinnode = c(10, 15, 20),
                           bag.fraction = c(.65, .8, 1),
                           optimal_trees = 0,          # a place to dump results
                           min_RMSE = 0                # a place to dump results
)                      

# total number of combinations
nrow(hyper_grid2)

#2nd grid search 
for(i in 1:nrow(hyper_grid2)) {
  # reproducibility
  set.seed(123)
  # train model
  gbm.tune2 <- gbm(
    formula = dgwl ~ .,
    distribution = "gaussian",
    data = random_gwla_train,
    n.trees = 5000,
    interaction.depth = hyper_grid2$interaction.depth[i],
    shrinkage = hyper_grid2$shrinkage[i],
    n.minobsinnode = hyper_grid2$n.minobsinnode[i],
    bag.fraction = hyper_grid2$bag.fraction[i],
    train.fraction = .75, #only uses 75% of the training data, the remaining 25% used for evaluation of performance
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )
  
  # add min training error and trees to grid
  hyper_grid2$optimal_trees[i] <- which.min(gbm.tune2$valid.error)
  hyper_grid2$min_RMSE[i] <- sqrt(min(gbm.tune2$valid.error))
}

hyper_grid2 %>% 
  dplyr::arrange(min_RMSE) %>%
  head(10)

#Final model
set.seed(123)
gbm.fit.final <- gbm(
  formula = dgwl ~ .,
  distribution = "gaussian",
  data = train_dataset,
  n.trees = 2461,
  interaction.depth = 10,
  shrinkage = 0.01,
  n.minobsinnode = 15,
  bag.fraction = 1.0, 
  train.fraction = 1,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
) 

print(gbm.fit.final)
sqrt(min(gbm.fit.final$valid.error))     #RMSE
sqrt(min(gbm.fit.final$train.error))



#visualize the model
##Variable importance
par(mar = c(5, 8, 1, 1))
summary(
  gbm.fit.final, 
  cBars = 10,
  method = relative.influence, # also can use permutation.test.gbm
  las = 2
)

##PDE plot
gbm.fit.final %>%
  partial(pred.var = "Discharge", n.trees = gbm.fit.final$n.trees, grid.resolution = 100) %>%
  autoplot(rug = TRUE, train = train_dataset) 

# get a few observations to perform local interpretation on Source:http://uc-r.github.io/lime
model_type.gbm <- function(x, ...) {
  return("regression")
}
predict_model.gbm <- function(x, newdata, ...) {
  pred <- predict(x, newdata, n.trees = x$n.trees)
  return(as.data.frame(pred))
}

local_obs <- test_dataset[1:2, ]

# apply LIME
explainer <- lime(train_dataset, gbm.fit.final)
explanation <- explain(local_obs, explainer, n_features = 5)
plot_features(explanation)

#Training metrics
pred <- predict(gbm.fit.final, n.trees = gbm.fit.final$n.trees, train_dataset)
rmse_train <- round(RMSE(pred, train_dataset$dgwl), 2)
RSquared_train <- round(R2_Score(pred, train_dataset$dgwl), 2)
plot(pred, train_dataset$dgwl, asp=1)
abline(coef=c(0,1), col="red", lwd=2)


#Testing
predTest <- predict(gbm.fit.final, n.trees = gbm.fit.final$n.trees, test_dataset)
rmse_tes <- round(RMSE(predTest, test_dataset$dgwl), 2)
RSquared_test <- round(R2_Score(predTest, test_dataset$dgwl), 2)
plot(predTest, test_dataset$dgwl, asp=1)
abline(coef=c(0,1), col="red", lwd=2)

RMSE <- c(rmse_train, rmse_tes)
Rsquared <- c(RSquared_train, RSquared_test)
Model <- c("Train", "Test")
metrics <- data.frame(Model, RMSE, Rsquared)
metrics


#Plot model simulated vs observed for a single well.
dgwlData2 <- read_csv("dgwlData.csv")
observedWell <- dgwlData2 %>% 
  filter(wellName=="Mt.Morris") 
observedWell <-observedWell[, -c(3, 16)] 
colnames(observedWell)[15] <- "dgwl"

#Imputed data
observed2 <- dgwlData2 %>% 
  filter(wellName=="Mt.Morris") %>% 
  imputate_na(DGWL_m, method = "mice") 

Simulated <- predict(gbm.fit.final, n.trees = gbm.fit.final$n.trees, observedWell)
plot(Simulated, observed2, asp=1)
abline(coef=c(0,1), col="red", lwd=2)
ggof(Simulated, Observed2)


unique(dgwlData2$wellName)

