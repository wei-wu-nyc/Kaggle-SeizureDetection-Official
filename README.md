# Kaggle-SeizureDetection-Official
My R code for Kaggle Competition Seizure Detection

I use the post (with a little editing) I posted to Kaggle group after the competition was over as this README. This post described the main ideas for my approach.

Congratulations to the winners. This competition turned out to be an interesting and fun competition for me unexpectedly.
I have learned and used many tips from the discussions from other Kagglers in this and other competitions. Even though I did not squeeze into top tier, I still would like to share some of my thoughts and approaches which may provide some aspects that have not been used/mentioned by others. The final private LB scores of this competition is so tight that if I had chose my best private score submission, (I chose my 2nd and 3rd best private scores, which is not too bad. The 2nd best is only 0.006 less than my best score.), I could have jumped by 4-5 positions. I hope someone can find some of the points useful. The sharing and learning environment at Kaggle is the best thing I like about Kaggle.

I did not plan to participate this competition as I saw the data set is pretty big seem like require a lot of time. Just about two weeks before the end, I saw that I could have time in my schedule to participate a competition and at around the same time it was announced that the daily submission limit was increased to 10 which make my late entrance more feasible. I thought I might take a crack at this.

I used R for this competition. The runs were mostly on my Home PC with i7 CPU and 12GB memory.

Overall consideration and approach:

a. Because the limited time available and the significant turn-around time to generate new set of features. I figured that I need to chose a set of features and stick to it. No time to do much feature selection/reduction.

b. I have no domain knowledge and did not participate in the first Seizue competition. Did not have much time to dive into the literatures. All the features I chose was from browsing through the forum of this and previous Seizure competition. And only choose items that I can understand and quickly implement.

c. One of the problem I see that I need to deal with at the beginning is the conflict of a performance metrics that includes all subjects and seperately trained models for each individual subjects. Seperately trained models for individual subjects may perform better for the individual subject, however, the probabilities generated from different models can have very different range, magnitude. And because for AUC metric, what matters is the ranking order of the predictions, when you combine results from individual models, because the differences in overall magnitudes or ranges of the predictions, simply combining the results from individual models may not yield best results. I decided at the beginning that for each type of model I use, I will have models trained for each subject using all the features available to that subject and a same type model that are global, trained using data from all subjects using a set of features common to all subjects. My plan was to later use the predictions from the global model to calibrate and combine the outputs of individual models. But I did not get time to explore this and at the end just simply ensemble both global and simply combined individual model results.

I will discuss my approach on the the following topics: Feature generation, Models used, Ensemble of different models.

Features:
a. To make features more comparable between subject and also to reduce the size, down sample signal data to 200Hz.

b. There were 3 evolving stages in the progress of feature generation in my approach. Each stage saw significant improvement in LB results.

c. Summary statistics only features: First, I saw in the post of the previous competition, some of the top winners only used summary statistics of the time series and FFT data. For each mat file, there are N number of channels. For each channel, besides the original signal series, I also generated delta and delta of delta series. For each of these series in time domain, I also generated FFT series. Then for each data series (time-domain and FFT) I generated a few basic summary statistics (mean, max, Stdev). Plus Frequency of FFT peak. 

Then also the overall summary statistics. First, generated statistic for a series of average of all channels. mean, max, stdev for all the quantities across all channels. Then the mean, max , SD and max,2nd max, 3rd max of elements of covariance matrix across all channel for both the time-domain series and FFT series.

Total 600+ features for dogs, little less number of common features, much more for patient_2. The best public LB results I could get with these features were: 0.657 for single model, 0.687 from ensemble

d. Summary statistics plus FFT: Adding FFT series is the natural next step. To reduce the size of FFT series, borrowed from the thinking of @ruai at "what-type-of-models-are-people-using", somewhat different parameters, using only the front portion (0.2) and smoothing the raw FFT to get only 24 points per FFT series. To arrive at the value of these parameters, I just simply plot the FFT series of a few files with different cutoff and averaging settings to get what I feel that are manageable and still kept enough details of the FFT series.
2000+ features for dogs, best single model: 0.69284, ensemble: 0.717
e. With the tips from this discussion "increase-number-of-samples-by-using-short-time-segments":

Break each of the original 10Min series into 10 non-overlap 60s series. The train data increased 10-fold. Used the same features as above for these shorter series. The eventual prediction for each original 10Min series was taken the average or certain percentile of the 10 predictions from the shorter series.
Best single model: 0.80514, Best Ensemble: 0.81803

Model used:
a. After the Higgs-Boson competition, xgboost become my favorite tool. Fast, effective.
b. For this data set SVM performed very well.
c. glmnet is competitive
d. also tried gbm and randomForest. These models had much worse results and very slow. I wish I did not waste too much time at the end to run these models instead I should have put more efforts into better tuning the better performed models.
My thought was to use as many types of fundamentally different models to ensemble as possible. In hindsight, I probably should ensemble several runs of the better performed models with different parameters. (Not tested.)
e. Each model type, both individual model and global model were constructed and tested.

Ensemble method:
a. AUC only depends on the order of the predictions. From previous competitions, I found the most effective way to ensemble a set of results from different models for AUC metric is to convert each solution to ranks and take the average ranks of all the solutions.
b. This ensemble method provides about 0.01-0.03 improvements.

Cross-Validation

2-fold CV with 10 shuffles, based on what the data description page used as "series", i.e. mat files with different sequence number but belong to the same series are in and out of folds together. The resulted AUC higher than public LB, but direction wise are pretty good.

My public and private LB scores are pretty much in sync. My public LB position was 25, private LB position jump to 11. I think my approach of constructing models trained both globally and individually helped in giving the stability of my results.

This is my experience. Hope you can find something useful.
