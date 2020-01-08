library(tidyverse)
library(VIM)
library(corrplot)
library(car)
library(EnvStats)
library(mlbench)
library(ggplot2)
library(tidyverse)
library(mice)
library(VIM)
library(forcats)
library(caret)
library(dplyr)
library(grid)
library(class)
library(cluster)
library(factoextra)
library(kohonen)
library(pROC)
library(ROCR)
library(rgl) 
#read csv file
Data <- read.csv(file="MIP3HWell.csv", header=TRUE)

#################################################################
#create Drilling and Geomechanical Histogram and Boxplots
#drillingdata reshaping dataframe
Data.df1 <- Data[,c(1,3,4,5,6,7,8,9,10,11,12,13,14)] %>% drop_na()
melt.df1 <- reshape2::melt(Data.df1,id='Hole.Depth')
Data.df2 <- Data[,c(1,15,16,17,18,19,20,21,22,23,24,25,26,27,28)] %>% drop_na()
melt.df2 <- reshape2::melt(Data.df2,id='Hole.Depth')

#Hitogram for drilling data
ggplot(data=melt.df1,aes(x=value,fill=variable)) +
  geom_histogram(bins=40)+facet_wrap(~variable,scale='free')+scale_x_continuous()
ggplot(data=melt.df2,aes(x=value,fill=variable)) +
  geom_histogram(bins=40)+facet_wrap(~variable,scale='free')

#Boxplots for drilling data
ggplot(data=melt.df1,aes(y=value,fill=variable)) +
  geom_boxplot()+facet_wrap(~variable,scale='free') 
ggplot(data=melt.df2,aes(y=value,fill=variable)) +
  geom_boxplot()+facet_wrap(~variable,scale='free') 

#geomechanical data reshaping dataframe
Data.gf1 <- Data[,c(1,29,30,31,34,35,36)] %>% drop_na()
melt.gf1 <- reshape2::melt(Data.gf1,id='Hole.Depth')

#Hitogram for geomechanical data
ggplot(data=melt.gf1,aes(x=value,fill=variable)) +
  geom_histogram(bins=40)+facet_wrap(~variable,scale='free')+scale_x_continuous()

#Boxplots for drilling data
ggplot(data=melt.gf1,aes(y=value,fill=variable)) +
  geom_boxplot()+facet_wrap(~variable,scale='free')+scale_x_continuous()

################################################################################################
#handling missingness in the data
Data %>% select_if(is.numeric) %>% mutate_all(is.na) %>% summarise_all(mean) #missingmess for each numeric attribute
par(mfrow=c(1,1))
md.pattern(Data, rotate.names = TRUE) #Show missing data pattern using MICE package
aggr(Data[,c(15:36)], cex.axis = 0.6, gap = 1, ylab = c("Proportion of Missingness", "Missingness Pattern"))

#drop NA value with mean more than 70%
Data <- Data %>% dplyr::select(-c(colnames(Data %>% mutate_all(is.na) %>% select_if(function(col) mean(col) > 0.7))))
md.pattern(Data, rotate.names = TRUE) #Show missing data pattern using MICE package
aggr(Data[,c(15:32)], cex.axis = 0.6, gap = 1, ylab = c("Proportion of Missingness", "Missingness Pattern"))
Data %>% mutate_all(is.na) %>% summarise_all(mean) 

#filling missingness value 
Data.imp <- mice(Data,m=10,maxit=20,method='norm',seed=100)
summary(Data.imp)
densityplot(Data.imp) #density plots show change is distribution of imputed data as compared to available data 
# but the distributions not match so we remove all missing values
#removing incomplete rows
Data = na.omit(Data)
Data %>% mutate_all(is.na) %>% summarise_all(mean) #all giving 0 mean of missingness
aggr(Data[,c(15:32)], cex.axis = 0.6, gap = 1, ylab = c("Proportion of Missingness", "Missingness Pattern"))

#########################################
#BUILD KMEANS CLUSTERS USING GEOMECH DATA 
#remove missing values
data_geomec <- Data %>% dplyr::select(c("Brittleness","YME","PR")) %>% drop_na()

#scaling/ standardizing data
geomec.scaled<-scale(data_geomec)

#determine optimal number of clusters 
  #quote by https://uc-r.github.io/kmeans_clustering
#determine optimal clusters basde on within-sum of squares and silhouette width
fviz_nbclust(data_geomec, kmeans, method = "wss")
fviz_nbclust(data_geomec, kmeans, method = "silhouette")

#create Kmeans classification
set.seed(40)
k1 <-kmeans(geomec.scaled,centers=3,nstart = 25)
fviz_cluster(k1, data =geomec.scaled )

storedcluster <- k1$cluster #stored the clusters
k1$centers

#Brittleness        YME         PR
#1  -0.1022563 -0.4771597  0.1449352
#2   4.5593713  3.3012232 -3.0077145
#3  -0.3215871  1.4480703 -0.1109248

#check classification with geomechanical parameters 
data_geomec$BrittlenessClassification <- as.factor(k1$cluster)
  #BOXPLOTS
ggplot(data_geomec,aes(x=BrittlenessClassification,y=YME,fill=BrittlenessClassification)) + geom_boxplot() + facet_wrap(~BrittlenessClassification,scale="free") +
    scale_y_continuous(limits = c(1,7.5))
ggplot(data_geomec,aes(x=BrittlenessClassification,y=PR,fill=BrittlenessClassification)) + geom_boxplot() + facet_wrap(~BrittlenessClassification,scale="free")+
  scale_y_continuous(limits = c(0,0.5))
  #PAIRPLOTS
ggplot(data_geomec,aes(x=PR,y=YME,color=BrittlenessClassification)) + geom_point() #+ facet_wrap(~BrittlenessClassification,scale="free")

#check classification with lithology 
data_litho <- Data %>% dplyr::select(c("GR","RHOB","DTCO")) %>% drop_na()
  #BOXPLOTS
data_litho$BrittlenessClassification <- as.factor(k1$cluster)
ggplot(data_litho,aes(x=BrittlenessClassification,y=GR,fill=BrittlenessClassification)) + geom_boxplot() + facet_wrap(~BrittlenessClassification,scale="free") +
  scale_y_continuous(limits = c(19,700))
ggplot(data_litho,aes(x=BrittlenessClassification,y=RHOB,fill=BrittlenessClassification)) + geom_boxplot() + facet_wrap(~BrittlenessClassification,scale="free")+
  scale_y_continuous(limits = c(2.4,2.9))
ggplot(data_litho,aes(x=BrittlenessClassification,y=DTCO,fill=BrittlenessClassification)) + geom_boxplot() + facet_wrap(~BrittlenessClassification,scale="free")+
  scale_y_continuous(limits = c(48,102))

#3D Visualization then PCA is applied, and then color code base on clusters
pc <- princomp(data_litho[,1:3], cor=TRUE, scores=TRUE)
summary(pc)

plot3d(pc$scores[,1:3], size=5, col=data_litho$BrittlenessClassification, main="k-means clusters")


#################################################
#cluster<- k1$cluster
#Data1  <- Data %>% drop_na(c("Brittleness","YME","PR","GR","RHOB","DTCO"))
Data$BrittlenessClassification <- as.factor(k1$cluster)

################################################

###############################################
  # Stratified random sample of data into train and test sets
train.index <- createDataPartition(Data$BrittlenessClassification,p=0.65,list=FALSE)#65% for train
train <- Data[train.index,]
test <- Data[-train.index,]

  #set out train data
geomec_train <- train %>% dplyr::select(c("Brittleness","YME","PR")) %>% drop_na()

################################################
#Unsupervised SOM
#the SOM training process.
#quote by: https://clarkdatalabs.github.io/soms/SOM_NBA
som_grid <- somgrid(xdim = 40, ydim=40, topo="hexagonal")
# Finally, train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(as.matrix(scale(geomec_train)), 
                 grid=som_grid, 
                 rlen=200, 
                 radius=2.5, 
                 keep.data = TRUE,
                dist.fcts= "euclidean")

# plot SOM
plot(som_model, main = "Default SOM Plot")
str(data_geomec) #recheck geomec dataframe

####################################
#Data prep for prediction (Classification)
#Normalization
normalize <- function(v){
  return((v-min(v))/(max(v)-min(v)))
}
#Create normalized dataset
BC = Data %>% dplyr::select(c(BrittlenessClassification))
Data_all = Data %>% dplyr::select(-c(Hole.Depth, Section, GR, RHOB, DTCO, Brittleness, YME, PR, BrittlenessClassification)) %>% glimpse
Data_all  = lapply(Data_all,normalize)

Data_all = cbind(Data_all, BC)

#data Splicing
set.seed(123)
trainIndex<- sample(1:nrow(Data_all),size=nrow(Data_all)*0.65,replace = FALSE) #random selection of 80% data.

train <- Data_all[trainIndex,] # 65% training data
train_BC = train %>% dplyr::select(BrittlenessClassification)
#train_knn = train_knn %>% dplyr::select(-c(BrittlenessClassification))

test <- Data_all[-trainIndex,] # remaining 35% test data
test_BC = test %>% dplyr::select(BrittlenessClassification)
test = test %>% dplyr::select(-c(BrittlenessClassification))

#######################################
#KNN
#build KNN model
set.seed(42)
control = trainControl(method = "repeatedcv", number = 10, repeats = 3, verboseIter = TRUE) 
knn_fit = train(BrittlenessClassification ~., data = train, 
                method = "knn", trControl = control, tuneLength = 20, preProcess = c("center", "scale"))
#k   Accuracy   Kappa    
#5  0.9820286  0.9481593
#7  0.9834691  0.9520598
#9  0.9842924  0.9543734
#11  0.9845669  0.9551447
#13  0.9844297  0.9547128
#15  0.9841554  0.9538519
#17  0.9838122  0.9527819
#19  0.9833319  0.9513921
#21  0.9832636  0.9511880
#23  0.9829208  0.9502380
#25  0.9827150  0.9496612
#27  0.9827150  0.9496802
#29  0.9827152  0.9496996
#31  0.9826466  0.9495411
#33  0.9828521  0.9501487
#35  0.9827150  0.9497314
#37  0.9831271  0.9509483
#39  0.9832643  0.9512926
#41  0.9836074  0.9522827
#43  0.9834702  0.9518383

#Accuracy was used to select the optimal model using the largest value.
#The final value used for the model was k = 11.

pred_knn = predict(knn_fit, newdata = test) #predict KNN into test set
caret::confusionMatrix(table(pred_knn, as.matrix(test_BC)))

#Confusion Matrix and Statistics

#pred_knn    1    2    3
#1 2037    2   37
#2    0   88    6
#3    6    0  442

#Accuracy : 0.9805          
#95% CI : (0.9745, 0.9855)
#No Information Rate : 0.7804          
#P-Value [Acc > NIR] : < 2.2e-16       
#Kappa : 0.9441          
#Mcnemar's Test P-Value : 1.165e-06

#ROC and AUC for prediction class 2
pred_knn_prob = predict(knn_fit, newdata = test, type = "prob")

plot(
  performance(
    prediction
    (pred_knn_prob[,2],factor(1*(test_BC == 2))),
    'tpr','fpr'
  )
)

auc_knn =   performance(
  prediction
  (pred_knn_prob[,2],factor(1*(test_BC == 2))),
  'auc'
)
auc_knn

#Slot "y.values":
#  [[1]]
#[1] 0.9940489

########################################
#GBM
set.seed(42)
control = trainControl(method = "repeatedcv", number = 3, repeats = 1, verboseIter = TRUE) 
hyper_grid = expand.grid(n.trees = c(200,400,800), interaction.depth = c(5,10,15), shrinkage = 0.01, n.minobsinnode = c(1,10))
gbm_fit = train(BrittlenessClassification ~., data = train, 
                method = "gbm", trControl = control, tuneGrid = hyper_grid, preProcess = c("center", "scale"))


hyper_grid_best = expand.grid(n.trees = c(800), interaction.depth = c(5), shrinkage = 0.01, n.minobsinnode = c(1))
#n.trees interaction.depth shrinkage n.minobsinnode
#1     800                 5      0.01              1
gbm_best = train(BrittlenessClassification ~., data = train, 
                 method = "gbm", trControl = control, tuneGrid = hyper_grid_best, preProcess = c("center", "scale")) 

#Stochastic Gradient Boosting 
#4860 samples
#24 predictor
#3 classes: '1', '2', '3' 
#Pre-processing: centered (24), scaled (24) 
#Resampling: Cross-Validated (3 fold, repeated 1 times) 
#Summary of sample sizes: 3241, 3240, 3239 
#Resampling results:
  
#Accuracy   Kappa    
#0.9862135  0.9601144

pred_gbm = predict(gbm_best, newdata = test)
pred_gbm_prob = predict(gbm_best, newdata = test, type = "prob")
caret::confusionMatrix(table(pred_gbm, as.matrix(test_BC)))

#Confusion Matrix and Statistics

#pred_gbm    1    2    3
#1 2037    2   32
#2    0   87    1
#3    6    1  452

#Overall Statistics

#Accuracy : 0.984           
#95% CI : (0.9784, 0.9884)
#No Information Rate : 0.7804          
#P-Value [Acc > NIR] : < 2.2e-16       
#Kappa : 0.954           
#Mcnemar's Test P-Value : 0.0001877 

#ROC and AUC for prediction class 2

plot(
  performance(
             prediction
                      (pred_gbm_prob[,2],factor(1*(test_BC == 2))),
             'tpr','fpr'
             )
     )

auc_gbm =   performance(
                  prediction
                       (pred_gbm_prob[,2],factor(1*(test_BC == 2))),
                  'auc'
                        )
auc_gbm
#Slot "y.values":
#  [[1]]
#[1] 0.9990594

############################################
#Neural Nets
#mlpML optimization and training
mlp_grid = expand.grid(layer1 = c(5,10,15,20),
                       layer2 = c(5,10,15,20),
                       layer3 = c(5,10,15,20))

mlp_fit = caret::train(BrittlenessClassification ~., data = train,
                       method = "mlpML", maxit = 100, learFunc = "Std_Backpropagation",learnFuncParams = c(0.01,0),
                       trControl = trainControl(method = "cv", number = 10,verboseIter = TRUE),
                       tuneGrid = mlp_grid, preProcess = c("center", "scale"),verbose = TRUE)

#best tuned model
mlp_grid_best = expand.grid(layer1 = 20,
                            layer2 = 20,
                            layer3 = 10)

#layer1 layer2 layer3
#1     20     20     10

mlp_best = caret::train(BrittlenessClassification ~., data = train,
                        method = "mlpML", maxit = 1000, learFunc = "Std_Backpropagation",learnFuncParams = c(0.01,0),
                        trControl = trainControl(method = "cv", number = 10,verboseIter = TRUE),
                        tuneGrid = mlp_grid_best, preProcess = c("center", "scale"),verbose = TRUE)

mlp_best

#Multi-Layer Perceptron, with multiple layers 

#4860 samples
#24 predictor
#3 classes: '1', '2', '3' 
#Pre-processing: centered (24), scaled (24) 
#Resampling: Cross-Validated (10 fold) 
#Summary of sample sizes: 4374, 4375, 4373, 4374, 4375, 4373, ... 
#Resampling results:
#Accuracy   Kappa    
#0.9864159  0.9604115
#Tuning parameter 'layer1' was held constant at a value of 20

pred_mlp = predict(mlp_best, newdata = test)
pred_mlp_prob = predict(mlp_best, newdata = test, type = "prob")
caret::confusionMatrix(table(pred_mlp, as.matrix(test_BC)))

#Confusion Matrix and Statistics

#pred_mlp    1    2    3
#1 2037    1   37
#2    3   89    1
#3    3    0  447

#Overall Statistics

#Accuracy : 0.9828           
#95% CI : (0.9771, 0.9874)
#No Information Rate : 0.7804          
#P-Value [Acc > NIR] : < 2.2e-16       
#Kappa : 0.9507           
#Mcnemar's Test P-Value : 8.923e-07 

#ROC and AUC for prediction class 2
plot(
  performance(
    prediction
    (pred_mlp_prob[,2],factor(1*(test_BC == 2))),
    'tpr','fpr'
  )
)

auc_gbm =   performance(
  prediction
  (pred_mlp_prob[,2],factor(1*(test_BC == 2))),
  'auc'
)

auc_gbm
#Slot "y.values":
#  [[1]]
#[1] 0.9953499



