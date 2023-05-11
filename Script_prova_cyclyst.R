#### Libraries ####

library(dplyr)
library(GEOquery)
library(genefilter)
library(randomForest)
library(caret) 
library(e1071) 
library(limma)
library(R.utils)
library(umap)
library(maptools)
library(MASS)
library(pROC)
setwd("C:\\Users\\filoa\\Desktop\\Network model")

set.seed(123)
#### Loading ####

cold <- getGEO("GSE85620")
cold_gene <- cold[[1]]
cold_data <- exprs(cold_gene)
anyNA(cold_data)

#### Naming ####

# Re-arranging names and labels, adding labels to order and then removing

part_cold <- cold[["GSE85620_series_matrix.txt.gz"]]@phenoData@data[["title"]]
Exposure <- as.vector(make.names(cold[["GSE85620_series_matrix.txt.gz"]]@phenoData@data[["group:ch1"]]))
Timing <- as.vector(cold[["GSE85620_series_matrix.txt.gz"]]@phenoData@data[["time:ch1"]])
cold_data <- as.data.frame(t(cold_data),
                        row.names = cold[["GSE85620_series_matrix.txt.gz"]]@phenoData@data[["title"]])
cold_data <- cold_data  %>%
  mutate(Timing = Timing, Exposure = Exposure) %>%
  arrange(Exposure, Timing)

Exposure_label <- cold_data$Exposure
Exposure <- cold_data$Exposure%in%"Cold.water.Immersion"
Placebo <- cold_data$Exposure%in%"Placebo"
Pre <- cold_data$Timing%in%"Pre"
Post <- cold_data$Timing%in%"Post"
cold_data <- cold_data[,1:19425]
View(cold_data[,19420:length(cold_data)])

#### PCA_1####

# PCA on all data set
pca <- cold_data[,1:19425] %>% prcomp

summary(pca) 
screeplot(pca)

paint_tot <- c(rep("skyblue",6), rep("magenta",6), rep("darkgreen",7),rep("darkorange",7))
paint_immertion <- c(rep("darkgreen",12),rep("darkorange",14))

plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=19, col = paint_tot)
legend("bottomleft", legend = c("Cold-Post", "Cold-Pre", "Ctrl-Post", "Ctrl-Pre"), col =unique(paint_tot), lwd = 4)

plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=19, col = paint_immertion)
legend("bottomleft", legend = c("Cold", "Ctrl"), col =unique(paint_immertion), lwd = 4)

#### Features selections ####

p_values_exposure <- c()
for (i in seq(1,19425)){
  t <- t.test(cold_data[Exposure,i],
              cold_data[Placebo,i])
  p_values_exposure <- append(p_values_exposure, t$p.value)
}
summary(p_values_exposure)
threshold <- 0.01
cold_meaningfull <- cold_data[,p_values_exposure < threshold] %>%  mutate(Exposure = Exposure_label)
hist(p_values_exposure)


tt40 <-rowttests(as.matrix(t(cold_data)), as.factor(Exposure_label))
ex3 <-cold_data[,which(tt40$p.value < threshold)]
hist(tt40$p.value)
cold <- cbind(as.data.frame(ex3),Exposure_label) 
colnames(cold)[ncol(cold)] <- "Exposure"

# confronting p_values
par(mfrow=c(1,2))
hist(p_values_exposure)
hist(tt40$p.value)
dev.off()

#### Predictions ####

n.controls <- 14
n.affected <- 12

train <- sample(1:(n.controls), (n.controls - 5)) 
test <- setdiff(1:(n.controls),train) 
test <- c(test, test + 12) 
train <-c(train, train + 12)

mod <- lda(Exposure ~ ., data = cold, prior = c(0.5,0.5), subset = train) 
plot(mod) 
mod.values <- predict(mod, cold[test,]) 
mod.values$class 
plot(mod.values$x[,1], ylab=c("LDA Axis"))
text(mod.values$x[,1], col = c(as.numeric(cold[train,"Exposure"]))) 

control <- trainControl(method = "cv", number = 10)
fit.lda <- train(Exposure~., data = cold, method = "lda", trControl = control, metric = metric) 
#varImp(fit.lda, scale = FALSE)

preds <- predict(mod, cold[test,])
preds$class 
table(preds$class, cold[test, "Exposure"])
metric <- "Accuracy"

roc_lda <- plot.roc(as.numeric(preds$class),as.numeric(as.factor(cold[test, "Exposure"])))
control <- trainControl(method = "repeatedcv", number = 10, repeats = 10) 
Ctrl <- trainControl(summaryFunction = twoClassSummary, classProbs = T)

fit.lda.2 <- train(Exposure ~ ., data = cold, method = "lda", metric = metric, trControl = control) 
fit.rf.2 <- train(Exposure ~ ., data = cold, method = "rf", metric = metric, trControl = control)
fit.m.lda2 <- train(Exposure ~ ., data = cold_meaningfull, method = "lda", metric = metric, trControl = control)
fit.m.rf2 <- train(Exposure ~ ., data = cold_meaningfull, method = "rf", metric = metric, trControl = control)

fit_LDA <- train(Exposure ~ ., data = cold, method = "lda",
                 metric = "ROC", preProc = c("center","scale"),
                 trControl = Ctrl)


results <- resamples(list(LDA = fit.lda.2, RF = fit.rf.2, M_lda = fit.m.lda2, M_rf = fit.m.rf2))
ggplot(results) + labs(y = "Accuracy")

varImp(fit_LDA)

#### Bho ####
fit_LDA <- lda(Exposure ~ ., data = cold)

#### PCA_2 ####

# after selecting the data
pca <- cold[,1:279] %>% prcomp

summary(pca) 
screeplot(pca)

plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=19, col = paint_tot)
legend("bottomleft", legend = c("Cold-Post", "Cold-Pre", "Ctrl-Post", "Ctrl-Pre"), col =unique(paint_tot), lwd = 4)

plot(pca$x[,1], pca$x[,2], xlab="PCA1", ylab="PCA2", main="PCA for components 1&2", type="p", pch=19, col = paint_immertion)
legend("topright", legend = c("Cold", "Ctrl"), col =unique(paint_immertion), lwd = 4)

#### Heat-maps ####

library(RColorBrewer)

app <- apply(cold[,1:279], 2, scale)
heatmap(app, col = terrain.colors(100))

#### Clustering ####

library(useful)

k <- 2
kmeans_result <- kmeans(cold[,1:279], k)
kmeans_result$totss
table(kmeans_result$cluster)

plot(kmeans_result, data = cold, main = "K-means Sub_meaning")+ 
  geom_text(aes(label =rownames(cold)), hjust = 0, vjust = 0)

dist_matrix <- dist(cold[,1:279])
hc_res <- hclust(dist_matrix, method = "ave")
k <- 2
groups <- cutree(hc_res, k = k)
table(groups)

plot(hc_res, hang = -1,labels = rownames(cold))
rect.hclust(hc_res, k = 2,which = NULL, x = NULL, h = NULL,cluster = NULL, border = 2)

#### Random Forest ####

## here we remove genes whose value is not > 0 in at least 20% of the samples 
ffun <- filterfun(pOverA(0.20,0.0))
t.fil <- genefilter(cold_data,ffun) 
small.eset <- log2(cold_data[t.fil,])
 
# small.eset<-log2(na.omit(e.mat)) 
dim(small.eset) 
y <- c(rep('C',12),rep('P',14))
rf <- randomForest(x = small.eset, as.factor(y), ntree=1000) 

# graph of sorted importance values 
plot(sort(rf$importance, decreasing=TRUE))    
# can also use: 
varImpPlot(rf) 
#extract the most 'important' genes 
probe.names <- rownames(rf$importance) 
top200 <- probe.names[order(rf$importance, decreasing=TRUE)[1:200]] 
# write.csv(top200, file = "probes-top200.txt", quote=FALSE, row.names= FALSE, col.names=FALSE)
plot(rf)

#### Glm-net ####

# TO CHECK AND TRY RIDGE

library(glmnet)
library(ROCR)

f <-  factor(y, labels = c("Cold_immertion","control"))
fit <- glmnet(cold, y, standardize=FALSE, family="binomial")

plot(fit, xvar= "lambda", label=TRUE) 
cfit <- cv.glmnet(as.matrix(cold[,1:279]),y,standardize=  FALSE, family="binomial") 

plot(cfit)
coef(cfit, s=cfit$lambda.min)

fit <- glmnet(cold[train,1:279], y[train], standardize = FALSE, family = "binomial") 
plot(fit) 

cfit <- cv.glmnet(as.matrix(cold[train,1:279]), y[train],standardize=FALSE, family="binomial")
plot(cfit)

predict(fit, as.matrix(cold[test,1:279]), type="class", s= cfit$lambda.min)

pred2 <-predict(fit, as.matrix(cold[test,1:279]), type="response", s=cfit$lambda.min) 
plot(performance(prediction(pred2, y[test]), 'tpr', 'fpr')) 


auc.tmp<-performance(prediction(pred2, y[test]), "auc") 
as.numeric(auc.tmp@y.values)

#### Lasso ####

library(caret)

control <-trainControl(method="cv", number=10)
metric <-"Accuracy"
fit.lasso<-train(Exposure ~., data = cold, method="glmnet",
                 family = "binomial",
                 tuneGrid= expand.grid(alpha = 1, lambda = seq(0,1,by=0.05)),
                 trControl= control,
                 metric = metric)
plot(fit.lasso) 

# comparison with other classification methods 
fit.lda <- train(Exposure~., data = cold, method="lda", metric=metric, trControl=control) 
fit.rf <- train(Exposure~., data = cold, method="rf", metric=metric, trControl=control) 
results <- resamples(list(RF=fit.rf, LDA=fit.lda, Lasso=fit.lasso)) 
summary(results) 
ggplot(results) + labs(y = "Accuracy")

#### Scudo ####

inTrain <- createDataPartition(f, list = FALSE) 
trainData <- t(cold[inTrain, 1:279]) 
testData <- t(cold[-inTrain, 1:279])
 
library(rScudo) 
trainRes <- scudoTrain(trainData, groups = f[inTrain],  nTop= 12, nBottom= 12, alpha = 0.05) 
trainRes 
 
upSignatures(trainRes)[1:5,1:5] 
consensusUpSignatures(trainRes)[1:5, ] 
 
trainNet <- scudoNetwork(trainRes, N = 0.2) 
scudoPlot(trainNet, vertex.label= NA) 
 
testRes <- scudoTest(trainRes, testData, f[-inTrain], nTop= 12, nBottom= 12) 
testNet <- scudoNetwork(testRes, N = 0.2) 
scudoPlot(testNet, vertex.label= NA) 
 
library(igraph) 

testClust <- igraph::cluster_spinglass(testNet, spins = 2) 
plot(testClust, testNet, vertex.label= NA) 
 
classRes<-scudoClassify(trainData, testData, N = 0.25,nTop= 12, nBottom= 12, trainGroups= f[inTrain], alpha = 0.5)
caret::confusionMatrix(classRes$predicted, f[-inTrain])

#### Gprofiling ####

# put the names of the gene selected in DAVID, paste the list, push gene list
# tell the identifier are in used
# options: background
# remove all the disease, leave the gene onthology, remove interactions

write.table(top200, paste("top200","selected",".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)



library(gprofiler2) 
# see vignette at https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html 
gostres<-gost(top200,
              organism = "hsapiens", ordered_query= FALSE, multi_query= FALSE,
              significant = TRUE, exclude_iea= FALSE, measure_underrepresentation= FALSE,
              evcodes= FALSE, user_threshold= 0.05, correction_method= "g_SCS",
              domain_scope= "annotated", custom_bg= NULL, numeric_ns= "",
              sources = NULL, as_short_link= FALSE)
names(gostres)
head(gostres$result)
# visualize results using a Manhattan plot
gostplot(gostres, capped = TRUE, interactive = TRUE)

# when ready,  create publication quality(static) plot + table of interesting terms/pathways 
p <-gostplot(gostres, capped = TRUE, interactive = FALSE) 
publish_gostplot(p, highlight_terms= c("GO:0048013", "REAC:R-HSA-3928663"),
                 width = NA, height = NA, filename = NULL)

#### selections variables ####

