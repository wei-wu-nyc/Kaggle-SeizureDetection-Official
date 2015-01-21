# This file is Wei Wu's code in participating
# Kaggle "American Epilepsy Society Seizure Prediction Challenge" Competition
# 
# Copyright (c) 2014, Wei Wu
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal 
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions: 
# 
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software. 
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
# IN THE SOFTWARE. 

library(e1071)
library(ROCR)

myauc <- function (pred, y) {
        myauc=as.numeric(performance(prediction(pred,y),"auc")@y.values)
}

get_periods <- function (a) {
	b=rep(0,length(a))
	for (i in 1:length(a)) {
		if (i == 1) {
			b[i]=1
		} else {
			if ((a[i]-a[i-1]) < 0) {
				b[i]=b[i-1]+1
			} else {
				b[i]=b[i-1]
			}
		}
	}
	return(b)
}

gen_k_folds <- function (k, pvals, nvals) {
	retval=list()
	pvals=sample(pvals, length(pvals), replace=FALSE)
	nvals=sample(nvals, length(nvals), replace=FALSE)
	for (i in 1:k) {
		retval[[i]]=list()
		retval[[i]][["pvals"]]=as.array(pvals[(round((i-1)*length(pvals)/k)+1):round(i*length(pvals)/k)])
		retval[[i]][["nvals"]]=as.array(nvals[(round((i-1)*length(nvals)/k)+1):round(i*length(nvals)/k)])
	}
	return(retval)
}

sub=c('Dog_1','Dog_2','Dog_3','Dog_4','Dog_5','Patient_1','Patient_2')
freq=200
dolog=1
addFFT=2
FFTratio=0.2
nsplit=10
holdout_types=c(1,2,3,4)    # c('Dog_1','Dog_2','Dog_3','Dog_4')

begTime=Sys.time()
alltrainmat=NULL
alltestmat=NULL
allholdoutmat=NULL
isub=1
allfields=list()
for (mysub in sub) {
	suppressWarnings(rm(trainmat, testmat, holdoutmat, summary_df))
	gc()
	DataPath=paste0('Data/',mysub)
	if (addFFT>=1) {
		freq_str=paste0(freq,"FFT",addFFT,"-",FFTratio)
	} else {
		freq_str=freq
	}
	if (nsplit > 1) {
		freq_str=paste0("Split",nsplit,"_Freq", freq_str)
	} else {
		freq_str=paste0("Freq",freq_str)
	}
	DataFile=paste0(DataPath,'/',mysub,freq_str,ifelse(dolog==1,"Log1",""),'.RData')
	cat(paste("Reading data from datafile:", DataFile,"\n")); flush.console()
	load(DataFile)
	nvar=ncol(trainmat)-4
	allfields[[isub]]=colnames(trainmat[,1:nvar])
	isub=isub+1
}
commonfields=allfields[[1]]
for (i in 1:length(allfields)) {
	commonfields=commonfields[commonfields %in% allfields[[i]]]
}
cat(paste("The number for common fields among all the subjects is:",length(commonfields),"\n"));flush.console()
oritraindata=list()
oritestdata=list()
oriholdoutdata=list()
isub=1
for (mysub in sub) {
	suppressWarnings(rm(trainmat, testmat, holdoutmat, summary_df))
	gc()
	DataPath=paste0('Data/',mysub)
	if (addFFT>=1) {
		freq_str=paste0(freq,"FFT",addFFT,"-",FFTratio)
	} else {
		freq_str=freq
	}
	if (nsplit > 1) {
		freq_str=paste0("Split",nsplit,"_Freq", freq_str)
	} else {
		freq_str=paste0("Freq",freq_str)
	}
	DataFile=paste0(DataPath,'/',mysub,freq_str,ifelse(dolog==1,"Log1",""),'.RData')
	cat(paste("Reading data from datafile:", DataFile,"\n")); flush.console()
	load(DataFile)
	summary_df=summary_df[grep(mysub,summary_df$fname),]
	extrafields=tail(colnames(trainmat),n=4)
	rownames(trainmat)=1:nrow(trainmat)
	rownames(testmat)=1:nrow(testmat)
        if (exists("holdoutmat")) {
	    rownames(holdoutmat)=1:nrow(holdoutmat)
        }
	trainmat=data.frame(trainmat)
	testmat=data.frame(testmat)
	trainmat$SubNum=isub
	testmat$SubNum=isub
        if (exists("holdoutmat")) {
	    holdoutmat=data.frame(holdoutmat)
	    holdoutmat$SubNum=isub
        }
#	trainmat=trainmat[,c(commonfields,"si", "SubNum",extrafields)]
#	testmat=testmat[,c(commonfields, "si", "SubNum")]
 	trainmat=trainmat[,c(commonfields, "SubNum",extrafields)]
 	testmat=testmat[,c(commonfields,  "SubNum", "si")]
        if (exists("holdoutmat")) {
	    holdoutmat=holdoutmat[,c(commonfields,  "SubNum", "si")]
        }
	nvar=ncol(trainmat)-5
	trainmat$period=isub*1000+get_periods(trainmat$seq)
	testmat$filename=rep(summary_df[grep('_test_',summary_df[,1]),1],each=nsplit)
	oritraindata[[isub]]=trainmat
	oritestdata[[isub]]=testmat
        if (exists("holdoutmat")) {
	    holdoutmat$filename=rep(summary_df[grep('_holdout_',summary_df[,1]),1],each=nsplit)
            oriholdoutdata[[isub]]=holdoutmat
        }
	if (is.null(alltrainmat)) {
		alltrainmat=trainmat
		alltestmat=testmat
                if (exists("holdoutmat")) {
		    allholdoutmat=holdoutmat
                }
	} else {
		alltrainmat=rbind(alltrainmat, trainmat)
		alltestmat=rbind(alltestmat, testmat)
                if (exists("holdoutmat")) {
		    allholdoutmat=rbind(allholdoutmat, holdoutmat)
                }
	}
	allfields[[isub]]=colnames(trainmat[,1:nvar])
	isub=isub+1
}
if (exists("holdoutmat")) {
    rm(holdoutmat)
}
rm(trainmat, testmat);gc()
pos_periods=unique(alltrainmat[alltrainmat$flag==1,'period'])
neg_periods=unique(alltrainmat[alltrainmat$flag==0,'period'])
alltrainmat$SubNum=as.factor(alltrainmat$SubNum)
alltestmat$SubNum=as.factor(alltestmat$SubNum)
if (!is.null(allholdoutmat)){
    allholdoutmat$SubNum=as.factor(allholdoutmat$SubNum)
}
cat(paste("Finished getting all the data files, Total time used:",format(Sys.time()-begTime),"\n"))
cat(paste("The sizes of the new train data and test data are:","(",paste(dim(alltrainmat),collapse=','),"),",
			"(",paste(dim(alltestmat),collapse=','),")\n"))
mysubnum=0
mycost=100
mygamma=0.001
nsplit=10
mydelflds=0

set.seed(20141109)
if (mysubnum == 1) {
	varnames=c(commonfields, "SubNum")
} else {
	varnames=commonfields
}
if (mydelflds == 1) {
	varnames=grep('^Del[12](|FFT)Cov',varnames,value=T,invert=T)
}


allsubtypes=as.numeric(unique(alltrainmat$SubNum))

inf_list=list()
all_pred_test=rep(NA,nrow(alltestmat))
if (!is.null(allholdoutmat)){
    all_pred_holdout=rep(NA,nrow(allholdoutmat))
}
for (mysubtype in allsubtypes) {
	begTime0=Sys.time()
	cat(paste("Training for Subject type",mysubtype,"...\n"));flush.console()
#	traindata=alltrainmat[alltrainmat$SubNum==mysubtype,]
	traindata=oritraindata[[mysubtype]]
	varnames=colnames(traindata)
	varnames=varnames[!varnames %in% c('seq','flag','period','ProbInd','SubNum','id','si.1')]
	if (mydelflds == 1) {
		varnames=grep('^Del[12](|FFT)Cov',varnames,value=T,invert=T)
	}
	svm.model=svm(traindata[,varnames], as.factor(traindata$flag),
		cost=mycost,
		gamma=mygamma,
#		epsilon=myepsilon, type=mytype, kernel=mykernel,
#		scale=ifelse(myscale==1, TRUE, FALSE),
		probability=TRUE)
#	testdata=alltestmat[alltestmat$SubNum==mysubtype,]
	testdata=oritestdata[[mysubtype]]
#	varnames=colnames(testdata)
#	varnames=varnames[!varnames %in% c('seq','flag','period','ProbInd','SubNum','id','si.1')]
	pred_test=predict(svm.model, testdata[,varnames], probability=TRUE)
	pred_test=attr(pred_test, "probabilities")[,2]
	all_pred_test[alltestmat$SubNum==mysubtype]=pred_test
        if (!is.null(allholdoutmat) && mysubtype %in% holdout_types) {
	    holdoutdata=oriholdoutdata[[mysubtype]]
	    pred_holdout=predict(svm.model, holdoutdata[,varnames], probability=TRUE)
	    pred_holdout=attr(pred_holdout, "probabilities")[,2]
	    all_pred_holdout[allholdoutmat$SubNum==mysubtype]=pred_holdout
        }
	if (length(pred_test[is.na(pred_test)]) > 0) {
	    stop("There is NA in the pred_test array.")
        }
	if ((!is.null(allholdoutmat)) && length(pred_holdout[is.na(pred_holdout)]) > 0) {
	    stop("There is NA in the pred_holdout array.")
        }

	endTime=Sys.time()
	cat(paste("Training for Subject type",mysubtype,"time used:", format(endTime-begTime0),"(",
			format(endTime-begTime),")\n"))
	flush.console()
}
if (length(all_pred_test[is.na(all_pred_test)]) > 0) {
	stop("There is NA in the prediction array.")
}
if ((!is.null(allholdoutmat)) && length(all_pred_holdout[is.na(all_pred_holdout)]) > 0) {
    stop("There is NA in the holdout prediction array.")
}

alltestmat$filenum=get_periods(alltestmat$si)
alltestmat$prob=all_pred_test
if (!is.null(allholdoutmat)) {
    allholdoutmat$filenum=get_periods(allholdoutmat$si)
    allholdoutmat$prob=all_pred_holdout
}

methods=c("mean","median","median")
percent_range=c(0,0.75,0.9)

for (i in 1:length(methods)) {
	avg_method=methods[i]
	percent=percent_range[i]
	if (avg_method == "mean") {
		tagg=aggregate(prob ~ filenum, data=alltestmat[,c('filenum','prob')], FUN=mean)
                if (!is.null(allholdoutmat)) {
		    tagg_holdout=aggregate(prob ~ filenum, data=allholdoutmat[,c('filenum','prob')], FUN=mean)
                }
	} else if (avg_method == "median") {
		tagg=aggregate(prob ~ filenum, data=alltestmat[,c('filenum','prob')], FUN=quantile, probs=c(percent))
                if (!is.null(allholdoutmat)) {
		    tagg_holdout=aggregate(prob ~ filenum, data=allholdoutmat[,c('filenum','prob')],
                                FUN=quantile, probs=c(percent))
                }
	} else {
		stop("Unknown averaging method.")
	}
	sub_df=data.frame(clip=alltestmat$filename[alltestmat$si==1],preictal=tagg$prob)
        if (!is.null(allholdoutmat)) {
	    holdout_df=data.frame(clip=allholdoutmat$filename[allholdoutmat$si==1],preictal=tagg_holdout$prob)
        }

	if (dolog==1) {
		freq_str0=paste0(freq_str,"Log1")
	}
	if (avg_method == "median") {
		avg_method=paste0(avg_method,percent)
	}
	filename=paste0("Official_SVM","_Individual_",avg_method,"_",freq_str0,"_cost",mycost,"_gamma",mygamma,
			"_delflds", mydelflds, "_sub.csv")
	write.csv(sub_df, file=filename, row.names=F, quote=F)
        if (!is.null(allholdoutmat)) {
	    filename=paste0("Holdout_SVM","_Individual_",avg_method,"_",freq_str0,"_cost",mycost,"_gamma",mygamma,
			"_delflds", mydelflds, "_sub.csv")
	    write.csv(holdout_df, file=filename, row.names=F, quote=F)
        }
}

cat(paste("Total time used:",format(Sys.time()-begTime),"\n"))
