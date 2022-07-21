# Library and data
library(rms)
library(pROC)
library(rmda)
train <-read.csv("E:/Experiments/YulanPeng/Wenwen/OSDFS/ScoreTrain-2OS.csv")
test <-read.csv("E:/Experiments/YulanPeng/Wenwen/OSDFS/ScoreTest.csv")

# Nomogram
dd=datadist(train)
options(datadist="dd")
f1 <- lrm(twoyOS~ twoOSCS
          + twoOSRS
          ,data = train,x = TRUE,y = TRUE)

nom <- nomogram(f1, fun=plogis,fun.at=c(.001, .01, seq(.1,.9, by=.4)), lp=F, funlabel="2-year Probability of OS")
plot(nom)

# ROC train
f2 <- glm(twoyOS~ twoOSCS
          + twoOSRS
          ,data = train,family = "binomial")

pre <- predict(f2, type='response')
plot.roc(train$twoyOS, pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, ci.type="bars", 
         of="thresholds",
         thresholds="best",
         print.thres="best",
         col="blue"
         #,identity=TRUE
         ,legacy.axes=TRUE,
         print.auc.x=ifelse(50,50),
         print.auc.y=ifelse(50,50)
)

# ROC test
pre1 <- predict(f2,newdata = test)
plot.roc(test$twoyOS, pre1,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, ci.type="bars",
         of="thresholds",
         thresholds="best",
         print.thres="best",
         col="blue",legacy.axes=TRUE,
         print.auc.x=ifelse(50,50),
         print.auc.y=ifelse(50,50)
)

# Calibration Curve train
rocplot1 <- roc(train$twoyOS, pre)
ci.auc(rocplot1)
cal <- calibrate(f1,  method = "boot", B = 1000)
plot(cal, xlab = "Nomogram Predicted Mutation", ylab = "Actual Mutation",main = "Calibration Curve")

# Calibration Curve test
rocplot2 <- roc(test$twoyOS,pre1)
ci.auc(rocplot2)
f3 <- lrm(test$twoyOS ~ pre1,x = TRUE,y = TRUE)
cal2 <- calibrate(f3,  method = "boot", B = 1000)
plot(cal2, xlab = "Nomogram Predicted Mutation", ylab = "Actual Mutation",main = "Calibration Curve")

# Decision Curve train
clinical_Rad<- decision_curve(twoyOS~ twoOSCS
                              + twoOSRS, data = train,
                              family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                              confidence.intervals= 0.95,study.design = 'case-control',
                              population.prevalence= 0.5)

List<- list(clinical_Rad)
plot_decision_curve(List,curve.names= c('Nomogram'),
                    cost.benefit.axis =FALSE,col = c('blue'),
                    confidence.intervals =FALSE,standardize = FALSE,
                    #legend.position = "none"
                    legend.position = "bottomleft"
)

# Decision Curve test
clinical_Rad1<- decision_curve(twoyOS~ twoOSCS
                               + twoOSRS, data = test,
                               family = binomial(link ='logit'), thresholds = seq(0,1, by = 0.01),
                               confidence.intervals= 0.95,study.design = 'case-control',
                               population.prevalence= 0.5)

List1<- list(clinical_Rad1)
plot_decision_curve(List1,curve.names= c('Nomogram'),
                    cost.benefit.axis =FALSE,col = c('blue'),
                    confidence.intervals =FALSE,standardize = FALSE,
                    legend.position = "bottomleft")
###
# Again.