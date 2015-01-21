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

require(plyr)

single_submissions=data.frame(subfile=c(
        "Official_XGBoost_median0.9_Split10_Freq200FFT2-0.2Log1_nround5000_eta0.05_depth6_lambda0.25_alpha0_child1_subsamp1_colsamp1_delflds1_sub.csv", # 0.80514 (Public LB)
        "Official_SVM_Individual_median0.9_Split10_Freq200FFT2-0.2Log1_cost100_gamma0.001_delflds0_sub.csv", # Old LB 0.79975,  New LB 0.79813
        "Official_SVM_median0.9_Split10_Freq200FFT2-0.2Log1_cost100_gamma0.001_delflds0_sub.csv", # Old LB 0.7879, New LB 0.78785
        "Official_glmnet_median0.75_Split10_Freq200FFT2-0.2Log1_alpha0_delflds0_sub.csv", # 0.70077
        "Official_glmnet_Individual_median0.9_Split10_Freq200FFT2-0.2Log1_alpha0_delflds0_sub.csv" #0.72466
	),
	subscore=c(0.80514, 0.79975, 0.78790, 0.70077, 0.72466),
	stringsAsFactors=FALSE)

### Setting do_holdout = 0 to get ensembled submission for test submission
### Setting do_holdout = 1 to get submission for post-competition holdout results
do_holdout=0

nmix=5
do_scale=2	# 0 no scaling, 1 - range scaling, 2 - rank scaling
weighted=5	# 0 - not weighted, 1 - score weighted, 2 - score^2, 3 - score^3
all_scores=NULL
weight_sum = 0
for (i in c(1:nmix)) {
	myfile=single_submissions$subfile[i]
        if (do_holdout == 1) {
            myfile=gsub("^Official_","Holdout_",myfile)
        }
	myscore=single_submissions$subscore[i]
	a=read.csv(myfile)
	if (do_scale == 1) {
		amax=max(a[,2])
		amin=min(a[,2])
		a[,2]=(a[,2]-amin)/(amax-amin)
	}
	if (do_scale == 2) {
#		a[order(a[,2]),2]=1:nrow(a)
		a[order(a[,2]),2]=(1:nrow(a))/nrow(a)
	}
	if (weighted == 1) {
		a[,2]=a[,2]*myscore
		weight_sum=weight_sum+myscore
	}
	if (weighted > 1) {
		a[,2]=a[,2]*myscore^weighted
		weight_sum=weight_sum+(myscore^weighted)
	}
	if (is.null(all_scores)) {
		all_scores=a
	} else {
		all_scores=join(all_scores, a, by="clip")
	}
}

if (weighted == 0) {
	sub_scores=data.frame(clip=all_scores$clip, preictal=rowMeans(all_scores[,2:ncol(all_scores)]))
} else {
	sub_scores=data.frame(clip=all_scores$clip, preictal=rowSums(all_scores[,2:ncol(all_scores)])/weight_sum)
}

if (do_holdout == 1) {
    fname=paste0("Holdoutmix_",nmix, "(xgb_svmInd_svm_glmnet_glmnetInd)", "_scale",do_scale,"_",
		ifelse(weighted,paste0("weighted",weighted),"equal"), "_sub.csv")
} else {
    fname=paste0("Officialmix_",nmix, "(xgb_svmInd_svm_glmnet_glmnetInd)", "_scale",do_scale,"_",
		ifelse(weighted,paste0("weighted",weighted),"equal"), "_sub.csv")
}
write.csv(sub_scores, row.names=F, quote=F, file=fname)
