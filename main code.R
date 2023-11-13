#Loads Libraries
library(pROC)
library(randomForest)
library(gtools)
library(vegan)
library(ade4)
library(randomForestExplainer)
library(effsize)
library(gtools)
library(RColorBrewer)

# iterative_rf_fixed_size -------------------------------------------------

iterative_rf_fixed_size = function(data,window1,window2,disease,control,iter,trainsize,testsize)
{
  set.seed(2023);
  featureProfile <- as.data.frame(matrix(NA,iter,ncol(data)));
  AUCArray <- NULL;
  SensitivityArray <- NULL;
  SpecificityArray <- NULL;
  numberDiseasedSamples <- NULL;
  numberControlSamples <- NULL;
  trainDiseaseSamples <- NULL;
  testDiseaseSamples <- NULL;
  trainControlSamples <- NULL;
  testControlSamples <- NULL;
  #threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
  for(i in 1:iter)
  {
    if(length(intersect(window1,disease)) > trainsize)
    {
      set.seed(i)
      trainDisease <- sample(intersect(window1,disease),trainsize,replace=FALSE);
    }
    else
    {
      set.seed(i)
      trainDisease <- sample(intersect(window1,disease),trainsize,replace=TRUE);
    }
    
    if(length(intersect(window1,control)) > trainsize)
    {
      set.seed(i)
      trainControl <- sample(intersect(window1,control),trainsize,replace=FALSE);
    }
    else
    {
      set.seed(i)
      trainControl <- sample(intersect(window1,control),trainsize,replace=TRUE);
    }
    
    tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
    rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
    TrainDiseaseTags <- NULL;
    TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
    TrainControlTags <- NULL;
    TrainControlTags[1:length(trainControl)] <- "Control";
    numberDiseasedSamples[i] <- nrow(trainDisease);
    numberControlSamples[i] <- nrow(trainControl);
    TrainTags <- c(TrainDiseaseTags,TrainControlTags);
    set.seed(i)
    rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
    #print("model created");
    featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
    if(length(setdiff(window1,window2)) == 0)
    {
      
      if(length(setdiff(intersect(window1,disease),trainDisease)) > testsize)
      {
        #print("LOOP1a");
        set.seed(i)
        testDisease <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=FALSE);
        #print(testDisease);
      }
      else
      {	
        #print(intersect(window1,disease));
        set.seed(i)
        testDisease <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=TRUE);
      }
      
      if(length(setdiff(intersect(window1,control),trainControl)) > testsize)
      {
        set.seed(i)
        testControl <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testControl <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=TRUE);
      }
      #print("test created");
      
    }
    else
    {
      if(length(intersect(window2,disease)) > testsize)
      {
        set.seed(i)
        testDisease <- sample(intersect(window2,disease),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testDisease <- sample(intersect(window2,disease),testsize,replace=TRUE);
      }
      
      if(length(intersect(window2,control)) > testsize)
      {
        set.seed(i)
        testControl <- sample(intersect(window2,control),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testControl <- sample(intersect(window2,control),testsize,replace=TRUE);
      }
      
    }
    tempTest <- rbind(data[testDisease,],data[testControl,]);
    TestDiseaseTags <- NULL;
    TestDiseaseTags[1:length(testDisease)] <- "Diseased";
    TestControlTags <- NULL;
    TestControlTags[1:length(testControl)] <- "Control";
    TestTags <- c(TestDiseaseTags,TestControlTags);
    rownames(tempTest) <- make.names(rownames(tempTest),unique=TRUE);
    rfTempPredict <- predict(rfTempComp,tempTest,type="vote",norm.votes=TRUE);
    print(i);
    AUCArray[i] <- auc(TestTags,rfTempPredict[,2])[1];
    print(median(AUCArray));
    SensitivityArray[i] <- length(which(predict(rfTempComp,tempTest[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTest[1:testsize,]));
    SpecificityArray[i] <- length(which(predict(rfTempComp,tempTest[(testsize+1):nrow(tempTest),])=="Control"))/length(predict(rfTempComp,tempTest[(testsize+1):nrow(tempTest),]));
    if(i == 100)
    {
      trainDiseaseSamples <- trainDisease;
      testDiseaseSamples <- testDisease;
      trainControlSamples <- trainControl;
      testControlSamples <- testControl;
    }
    i <- i + 1;
  }
  colnames(featureProfile) <- colnames(tempTrain);
  #returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
  returnList <- list("AUC"=AUCArray,"Sensitivity"=SensitivityArray,"Specificity"=SpecificityArray,"numberDiseasedSamples"=numberDiseasedSamples,"numberControlSamples"=numberControlSamples,"featureProfile"=featureProfile,"trainDisease100"=trainDiseaseSamples,"trainControl100"=trainControlSamples,"testDisease100"=testDiseaseSamples,"testControl100"=testControlSamples);
  return(returnList);
  
}


# iterative_rf_multiple_with_null_distribution ----------------------------

iterative_rf_multiple_with_null_distribution = function(data,window1,window2,disease,control,training_iter,testing_iter,trainsize,testsize)
{
  set.seed(2023);
  featureProfile <- as.data.frame(matrix(NA,training_iter,ncol(data)));
  AUCArraySelf <- matrix(NA,training_iter,testing_iter)
  SensitivityArraySelf <- matrix(NA,training_iter,testing_iter)
  SpecificityArraySelf <- matrix(NA,training_iter,testing_iter)
  AUCArrayNonSelf <- matrix(NA,training_iter,testing_iter)
  SensitivityArrayNonSelf <- matrix(NA,training_iter,testing_iter)
  SpecificityArrayNonSelf <- matrix(NA,training_iter,testing_iter)
  AUCArrayNull1 <- matrix(NA,training_iter,testing_iter)
  SensitivityNull1 <- matrix(NA,training_iter,testing_iter)
  SpecificityNull1 <- matrix(NA,training_iter,testing_iter)
  AUCArrayNull1 <- matrix(NA,training_iter,testing_iter)
  SensitivityArrayNull1 <- matrix(NA,training_iter,testing_iter)
  SpecificityArrayNull1 <- matrix(NA,training_iter,testing_iter)
  AUCArrayNull2 <- matrix(NA,training_iter,testing_iter)
  SensitivityArrayNull2 <- matrix(NA,training_iter,testing_iter)
  SpecificityArrayNull2 <- matrix(NA,training_iter,testing_iter)
  #threshold <- ifelse(length(intersect(window,disease)) <= 20,10,20);
  for(i in 1:training_iter)
  {
    if(length(intersect(window1,disease)) > trainsize)
    {
      set.seed(i)
      trainDisease <- sample(intersect(window1,disease),trainsize,replace=FALSE);
    }
    else
    {
      set.seed(i)
      trainDisease <- sample(intersect(window1,disease),trainsize,replace=TRUE);
    }
    
    if(length(intersect(window1,control)) > trainsize)
    {
      set.seed(i)
      trainControl <- sample(intersect(window1,control),trainsize,replace=FALSE);
    }
    else
    {
      set.seed(i)
      trainControl <- sample(intersect(window1,control),trainsize,replace=TRUE);
    }
    tempTrain <- rbind(data[trainDisease,],data[trainControl,]);
    rownames(tempTrain) <- make.names(rownames(tempTrain),unique=TRUE);
    TrainDiseaseTags <- NULL;
    TrainDiseaseTags[1:length(trainDisease)] <- "Diseased";
    TrainControlTags <- NULL;
    TrainControlTags[1:length(trainControl)] <- "Control";
    TrainTags <- c(TrainDiseaseTags,TrainControlTags);
    set.seed(i)
    rfTempComp <- randomForest(as.factor(TrainTags)~.,tempTrain);
    #print("model created");
    featureProfile[i,] <- sapply(colnames(tempTrain),function(x)(ifelse(x %in% rownames(rfTempComp$importance),rfTempComp$importance[x,],0)));
    for(k in 1:testing_iter)
    {
      if(length(setdiff(intersect(window1,disease),trainDisease)) > testsize)
      {
        set.seed(i)
        testDiseaseSelf <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testDiseaseSelf <- sample(setdiff(intersect(window1,disease),trainDisease),testsize,replace=TRUE);
      }
      
      if(length(setdiff(intersect(window1,control),trainControl)) > testsize)
      {
        set.seed(i)
        testControlSelf <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testControlSelf <- sample(setdiff(intersect(window1,control),trainControl),testsize,replace=TRUE);
      }
      
      if(length(intersect(window2,disease)) > testsize)
      {
        set.seed(i)
        testDiseaseNonSelf <- sample(intersect(window2,disease),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testDiseaseNonSelf <- sample(intersect(window2,disease),testsize,replace=TRUE);
      }
      
      if(length(intersect(window2,control)) > testsize)
      {
        set.seed(i)
        testControlNonSelf <- sample(intersect(window2,control),testsize,replace=FALSE);
      }
      else
      {
        set.seed(i)
        testControlNonSelf <- sample(intersect(window2,control),testsize,replace=TRUE);
      }
      tempTestSelf <- rbind(data[testDiseaseSelf,],data[testControlSelf,]);
      TestDiseaseTags <- NULL;
      TestDiseaseTags[1:length(testDiseaseSelf)] <- "Diseased";
      TestControlTags <- NULL;
      TestControlTags[1:length(testControlSelf)] <- "Control";
      TestTags <- c(TestDiseaseTags,TestControlTags);
      rownames(tempTestSelf) <- make.names(rownames(tempTestSelf),unique=TRUE);
      rfTempPredictSelf <- predict(rfTempComp,tempTestSelf,type="vote",norm.votes=TRUE);
      print(i)
      print(k)
      AUCArraySelf[i,k] <- auc(TestTags,rfTempPredictSelf[,2])[1];
      print(AUCArraySelf[i,k])
      SensitivityArraySelf[i,k] <- length(which(predict(rfTempComp,tempTestSelf[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTestSelf[1:testsize,]));
      SpecificityArraySelf[i,k] <- length(which(predict(rfTempComp,tempTestSelf[(testsize+1):nrow(tempTestSelf),])=="Control"))/length(predict(rfTempComp,tempTestSelf[(testsize+1):nrow(tempTestSelf),]));
      
      tempTestNonSelf <- rbind(data[testDiseaseNonSelf,],data[testControlNonSelf,]);
      #print(testDiseaseSelf)
      TestDiseaseTags <- NULL;
      TestDiseaseTags[1:length(testDiseaseNonSelf)] <- "Diseased";
      TestControlTags <- NULL;
      TestControlTags[1:length(testControlNonSelf)] <- "Control";
      TestTags <- c(TestDiseaseTags,TestControlTags);
      rownames(tempTestNonSelf) <- make.names(rownames(tempTestNonSelf),unique=TRUE);
      rfTempPredictNonSelf <- predict(rfTempComp,tempTestNonSelf,type="vote",norm.votes=TRUE);
      #print(k);
      AUCArrayNonSelf[i,k] <- auc(TestTags,rfTempPredictNonSelf[,2])[1];
      print(AUCArrayNonSelf[i,k])
      SensitivityArrayNonSelf[i,k] <- length(which(predict(rfTempComp,tempTestNonSelf[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTestNonSelf[1:testsize,]));
      SpecificityArrayNonSelf[i,k] <- length(which(predict(rfTempComp,tempTestNonSelf[(testsize+1):nrow(tempTestNonSelf),])=="Control"))/length(predict(rfTempComp,tempTestNonSelf[(testsize+1):nrow(tempTestNonSelf),]));
      set.seed(i)
      testDiseaseNull1 <- sample(sample(c(testDiseaseSelf,testDiseaseNonSelf)),testsize,replace=FALSE)
      set.seed(i)
      testControlNull1 <- sample(sample(c(testControlSelf,testControlNonSelf)),testsize,replace=FALSE)
      
      tempTestNull1 <- rbind(data[testDiseaseNull1,],data[testControlNull1,]);
      TestDiseaseTags <- NULL;
      TestDiseaseTags[1:length(testDiseaseNull1)] <- "Diseased";
      TestControlTags <- NULL;
      TestControlTags[1:length(testControlNull1)] <- "Control";
      TestTags <- c(TestDiseaseTags,TestControlTags);
      rownames(tempTestNull1) <- make.names(rownames(tempTestNull1),unique=TRUE);
      rfTempPredictNull1 <- predict(rfTempComp,tempTestNull1,type="vote",norm.votes=TRUE);
      #print(k);
      AUCArrayNull1[i,k] <- auc(TestTags,rfTempPredictNull1[,2])[1];
      print(AUCArrayNull1[i,k])
      SensitivityArrayNull1[i,k] <- length(which(predict(rfTempComp,tempTestNull1[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTestNull1[1:testsize,]));
      SpecificityArrayNull1[i,k] <- length(which(predict(rfTempComp,tempTestNull1[(testsize+1):nrow(tempTestNull1),])=="Control"))/length(predict(rfTempComp,tempTestNull1[(testsize+1):nrow(tempTestNull1),]));
      set.seed(i+1)
      testDiseaseNull2 <- sample(sample(c(testDiseaseSelf,testDiseaseNonSelf)),testsize,replace=FALSE)	
      set.seed(i+1)
      testControlNull2 <- sample(sample(c(testControlSelf,testControlNonSelf)),testsize,replace=FALSE)
      
      tempTestNull2 <- rbind(data[testDiseaseNull2,],data[testControlNull2,]);
      TestDiseaseTags <- NULL;
      TestDiseaseTags[1:length(testDiseaseNull2)] <- "Diseased";
      TestControlTags <- NULL;
      TestControlTags[1:length(testControlNull2)] <- "Control";
      TestTags <- c(TestDiseaseTags,TestControlTags);
      rownames(tempTestNull2) <- make.names(rownames(tempTestNull2),unique=TRUE);
      rfTempPredictNull2 <- predict(rfTempComp,tempTestNull2,type="vote",norm.votes=TRUE);
      #print(k);
      AUCArrayNull2[i,k] <- auc(TestTags,rfTempPredictNull2[,2])[1];
      SensitivityArrayNull2[i,k] <- length(which(predict(rfTempComp,tempTestNull2[1:testsize,])=="Diseased"))/length(predict(rfTempComp,tempTestNull2[1:testsize,]));
      SpecificityArrayNull2[i,k] <- length(which(predict(rfTempComp,tempTestNull2[(testsize+1):nrow(tempTestNull2),])=="Control"))/length(predict(rfTempComp,tempTestNull2[(testsize+1):nrow(tempTestNull2),]));	
      
      k <- k + 1;
    }	
    
    i <- i + 1;
  }
  colnames(featureProfile) <- colnames(tempTrain);
  #returnList <- list("featureProfile"=featureProfile,"AUC"=AUCArray,"Accuracy"=AccuracyArray);
  returnList <- list("AUCSameAge"=AUCArraySelf,"SensitivitySameAge"=SensitivityArraySelf,"SpecificitySameAge"=SpecificityArraySelf,"AUCDiffAge"=AUCArrayNonSelf,"SensitivityDiffAge"=SensitivityArrayNonSelf,"SpecificityDiffAge"=SpecificityArrayNonSelf,"AUCNull1"=AUCArrayNull1,"SensitivityNull1"=SensitivityArrayNull1,"SpecificityNull1"=SpecificityArrayNull1,"AUCNull2"=AUCArrayNull2,"SensitivityNull2"=SensitivityArrayNull2,"SpecificityNull2"=SpecificityArrayNull2,"featureProfile"=featureProfile);
  return(returnList);
  
}

# rank_scale --------------------------------------------------------------

rank_scale=function(x)
{
  x <- rank(x);
  y <- (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)));
  return(y);
}

# wilcox_batch ------------------------------------------------------------

wilcox_batch = function(x,y)
{
  p_array <- NULL;
  type_array <- NULL;
  mean1_array <- NULL;
  mean2_array <- NULL;
  x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
  y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
  z <- intersect(rownames(x),rownames(y));
  for(i in 1:length(z))
  {
    p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
    type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
    mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
    mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
    i <- i + 1;
  }
  out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array),mean1_array,mean2_array));
  rownames(out) <- z;
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  return(t(out));
}

# batch_dunns -------------------------------------------------------------

batch_dunns <- function(data, groups) {
  num_batches <- dim(data)[2]
  kruskal_p <- NA
  name_all <- colnames(data)
  for (i in 1:num_batches) {
    # Perform Kruskal-Wallis test
    kruskal_result <- kruskal.test(data[,name_all[i]], groups)
    # Extract p-values
    kruskal_p[i] <- kruskal_result$p.value
  }
  CorrectedP <- p.adjust(kruskal_p,method = "BH")
  # Return the results as a list
  out <- as.data.frame(cbind(kruskal_p,CorrectedP))
  rownames(out) <- colnames(data);
  out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
  return(data.frame(t(out)))
}














