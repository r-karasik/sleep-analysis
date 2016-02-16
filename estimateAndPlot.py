import matplotlib.pyplot as plt
import matplotlib
import numpy
import scipy.stats
from sklearn.linear_model import LogisticRegression
from sklearn import datasets, linear_model
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn import cross_validation
from sklearn import cluster
from sklearn.svm import SVR
import sklearn
from sklearn.manifold import MDS
import scipy.spatial
import random
import psycopg2

import parse
import analysis
from parse import dateDiff


#**********************************
#*** estimation
def dataForInitialAnalysis(children, analysis):
    #generates data for machine learning with "initial" data (prior history on sleep disorder) provided from the users
    index=analysis.keys()
    random.shuffle(index)
    #print "Does shuffle work", index
    #collect data in machine learning format
    X=[]
    Y=[]
    for i in index:
        init=children[i].initial
        if analysis[i].start_date:
            s_date=dateDiff(start, analysis[i].start_date)
            #print i, s_date, analysis[i].start_date
            X.append([init.month_count, init.night_count, s_date])
            Y.append(len(analysis[i].use))
    X=numpy.array(X)
    Y=numpy.array(Y)
    print 'mean', numpy.mean(X[:,0])
    print 'sdt', numpy.std(X[:,0])
    print 'sample X', X[:10]
    plotInitData(X,Y)
    return index, X, Y

def plotInitData(X,Y):
    #plots initial data to see if there are any obvious correlations
    plt.figure(1)
    plt.title('per month count vs len use')
    plt.xlabel('per month count')
    plt.ylabel('len use')
    plt.plot(X[:,0], Y, '*')
    #plt.savefig('figures/result/monthCountVsUse.jpg')
    plt.figure(2)
    plt.title('per night count vs len use')
    plt.xlabel('per night count')
    plt.ylabel('len use')
    plt.plot(X[:,1], Y, '*')
    #plt.savefig('figures/result/nightCountVsUse.jpg')
    plt.figure(3)
    plt.title('start date vs len use')
    plt.xlabel('start date')
    plt.ylabel('len use')
    plt.plot(X[:,2], Y, '*')
    #plt.savefig('figures/result/startDvsLenUse.jpg')
    plt.show()


def estLenUseFromInit(children, analysis, start):
    #estimates length of use from initial conditions
    lr= LogisticRegression(C=0.5)
    index, X, Y=dataForInitialAnalysis(children, analysis)
    total=len(index)
    train_num=int(total*.5) #number of data points used for train  set
    cr_val_num=int(total*.25) #number of points used for cross correlation set
    #print 'partition stats:', total, train_num, cr_val_num
    #normalize data for machine learning
    X=sklearn.preprocessing.normalize(1.0* X)
    print 'sample X after renom', X[:10]
    #print numpy.shape(X)
    #print 'remove last component', numpy.shape(X[:, :2])
    X_train=X[:train_num]
    X_cr_val=X[train_num:train_num+cr_val_num]
    X_test=X[train_num+cr_val_num:]
    Y_train=Y[:train_num]
    Y_cr_val=Y[train_num:train_num+cr_val_num]
    Y_test=Y[train_num+cr_val_num:]
    Y_tr=dict()
    Y_crv=dict()
    Y_pred_tr=dict()
    Y_pred_cv=dict()
    score=[]
    J_cv=[]
    J_tr=[]
    print 'X train char', numpy.shape(X_train)
    maxD=20
    for k in range(1,maxD):
        Y_tr[k]=Y_train>k
        Y_crv[k]=Y_cr_val>k
        
        lr.fit(X_train, Y_tr[k])
        print 'd=', k
        print lr.coef_
        print 'score', lr.score(X_cr_val, Y_crv[k]), lr.score(X_train, Y_tr[k])
        Y_pred_tr=lr.predict(X_train)
        Y_pred_cv=lr.predict(X_cr_val)
        print 'true', numpy.sum(Y_tr[k]), numpy.sum(Y_crv[k])
        print numpy.sum(Y_pred_tr), numpy.sum(Y_pred_cv)
        y_tr=1*(Y_pred_tr-Y_tr[k])
        J_tr.append((y_tr.dot(y_tr))/(2.*len(y_tr)))
        y_cv=1*(Y_pred_cv-Y_crv[k])
        #print y_cv, y_cv.dot(y_cv)
        J_cv.append((y_cv.dot(y_cv))/(2.*len(y_cv)))
        score.append(lr.score(X_cr_val, Y_crv[k]))
        #print 'probability', lr.predict_proba(X_cr_val)
    #print 'error', J_cv, J_tr
    plt.figure()
    plt.title('cross validation scores for different num of use days')
    plt.ylabel('score')
    plt.xlabel('number of days app was used')
    plt.plot(range(1,maxD), J_tr, 'b')
    plt.plot(range(1,maxD), J_cv, 'g')
    #plt.savefig('figures/result/initCondPred.jpg')
    plt.show()

######***********************************
#******Machine Learning analysis for long-term users
##****************************************

def dataForML(analysis):
    #prepares data to be used for machine learning analysis
    index=analysis.keys()
    random.shuffle(index) #shuffle users to ensure that data for ML is fair (no hidden trends based on when users joined the system)
    #print "Does shuffle work", index

    #initialize machine learning variables
    X=[]
    Y=[]
    Y2=[]
    Y3=[]
    Y4=[]
    #print "problem with analysis data structure?"
    for cid in index:
        ch=analysis[cid]
        #print cid, ch.std_time_for_bed

        if len(ch.biggest_improvement)>0:
            vector=[ch.std_time_for_bed, ch.std_event1_time, 
                    ch.std_time_to_event1, ch.avg_time_to_event1, 
                    ch.std_offset_for_device_start, ch.avg_offset_for_device_start, 
                    ch.std_device_start_time, ch.std_device_interval, 
                    ch.use_rate, len(ch.use), 
                    len(ch.vibration_use)]
            #vector=[ch.use_rate]
            #print 'time_for_bed?', cid, ch.std_time_for_bed
            
            valid=True
            for el in vector:
                if numpy.isnan(el):
                    #print 'child with id=', cid, 'invalid data point is', vector
                    valid=False

            if ch.corr_dev_reactgood is None or numpy.isnan(ch.corr_dev_reactgood):
                valid = False

            if valid:
                X.append(vector)
                Y.append(ch.biggest_improvement[2])
                Y2.append(ch.biggest_improvement[0])
                Y3.append(ch.corr_dev_reactgood)
                Y4.append(int(ch.works_flag))


    return X, Y, Y2, Y3, Y4


def plotMDS(X, Y):
    #computes and plots MDS (measure for how well data separates)
    D = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(X))
    tmodel = MDS(n_components=2, dissimilarity='precomputed')
    X2D = tmodel.fit_transform(D)
    plt.figure()
    plt.title('MDS')
    plt.ylabel('MDS1')
    plt.xlabel('MDS2')
    plt.scatter(X2D[:, 0], X2D[:, 1], c=Y)
    plt.show()

def plotPCA(X, Y, comment):
    #computes and plots PCA (another measure for how well data separates)
    X_Y = numpy.hstack((X, numpy.transpose([Y])))
    print X.shape, Y.shape, X_Y.shape
    for nc in range(1, 12):
        pca = PCA(n_components=nc, whiten=True)
        trans = pca.fit_transform(X)
        print "trying PCA with nc", nc, ", ratio is", numpy.sum(pca.explained_variance_ratio_)
        print pca.explained_variance_ratio_
        if nc == 2:
            print trans.shape, Y.shape
            plt.scatter(trans[:, 0], trans[:, 1], c=Y, edgecolor='face')
            plt.colorbar()
            plt.title(comment)
            plt.xlabel('PCA1')
            plt.ylabel('PCA2')
    plt.show()



def randomForestRegressorStudy(X,Y, setSize, comment):
    #runs random forest regressor on the data to see the performance of the prediction and to determine predictive features 
    X_train=X[:setSize]
    X_test=X[setSize:]
    Y_train=Y[:setSize]
    Y_test=Y[setSize:]

    rf_reg=RandomForestRegressor(n_estimators=10)
    rf_reg.fit(X_train, Y_train)
    Y_pred=rf_reg.predict(X_train)
    print "random forest regressor for "+comment, rf_reg.score(X_train, Y_train), rf_reg.score(X_test, Y_test)
    print "feature importances", rf_reg.feature_importances_

    scores = cross_validation.cross_val_score(rf_reg, X, Y, cv=5)
    print "cross-validation"
    print scores

def linearRegressionStudy(X, Y, setSize, comment):
    #runs linear regression on the data
    X_train=X[:setSize]
    X_test=X[setSize:]
    Y_train=Y[:setSize]
    Y_test=Y[setSize:]

    regr = linear_model.LinearRegression()
    regr.fit(X_train, Y_train)
    Y_pred_train = regr.predict(X_train)
    Y_pred_test = regr.predict(X_test)

    print "linear regression for "+comment
    # The coefficients
    print('Coefficients: %s' % regr.coef_)
    # The mean square error
    print("Residual sum of squares: %.2f"
          % numpy.mean((Y_pred_test - Y_test) ** 2))
    # Explained variance score: 1 is perfect prediction
    print('Variance score: %.2f' % regr.score(
        X_test, Y_test))
    print "train score", regr.score(X_train, Y_train)

    scores = cross_validation.cross_val_score(regr, X, Y, cv=10)
    print "cross-validation"
    print scores


def svmRegressorStudy(X, Y, setSize, comment):
    #runs svm regressor on the data
    X_train=X[:setSize]
    X_test=X[setSize:]
    Y_train=Y[:setSize]
    Y_test=Y[setSize:]

    svm=SVR()
    svm.fit(X_train, Y_train)
    print 'svm regressor '+comment
    s1 = svm.score(X_train, Y_train)
    s2 = svm.score(X_test, Y_test)
    print 'svm score for ', s1, s2

def logRegressionStudy(X, Y, setSize, comment):
    #runs logistic regression (classification) on the data
    X_train=X[:setSize]
    X_test=X[setSize:]
    Y_train=Y[:setSize]
    Y_test=Y[setSize:]


    score=[]
    J_test=[]
    J_tr=[]
      
    lr= LogisticRegression(C=0.5)  
    lr.fit(X_train, Y_train)
    print 'logistic regression '+comment, lr.score(X_train, Y_train), lr.score(X_test, Y_test)
    print 'coeff', lr.coef_

    Y_pred_train=lr.predict(X_train)
    Y_pred_test=lr.predict(X_test)
    #print Y_pred_tr, Y_train
    #print numpy.shape(Y_train), numpy.shape(Y_pred_tr)

    y_tr=1*(Y_pred_train-Y_train)
    J_tr.append((y_tr.dot(y_tr))/(2.*len(y_tr)))
    y_test=1*(Y_pred_test-Y_test)
    #print y_cv, y_cv.dot(y_cv)
    J_test.append((y_test.dot(y_test))/(2.*len(y_test)))
    score.append(lr.score(X_test, Y_test))
    print 'logistic regression '+comment, J_tr, J_test, score
    
    scores = cross_validation.cross_val_score(lr, X, Y, cv=8)
    print "cross-validation"
    print scores


def randomForestClassifierStudy(X, Y, setSize, comment):
    #runs random forrest classifier (classification) on the data
    X_train=X[:setSize]
    X_test=X[setSize:]
    Y_train=Y[:setSize]
    Y_test=Y[setSize:]


    cols = range(X_train.shape[1])
    print 'random forest with disappearing data'
    cvdata = []
    cvdata_avg = []
    cvlabels = []
    while len(cols) > 1:
        lr=RandomForestClassifier(n_estimators=5,
                                  class_weight='balanced')
        npcols = numpy.array(cols)
        lr.fit(X_train[:, npcols], Y_train)            
        coef = lr.feature_importances_
        maxi = numpy.argmax(coef)
        print 'cols=%d drop=%d train=%.4f test=%.3f %s' % (
            len(cols), cols[maxi], 
            lr.score(X_train[:, npcols], Y_train), 
            lr.score(X_test[:, npcols], Y_test),
            ','.join('%.4f' % r for r in coef))
            
        scores = cross_validation.cross_val_score(lr, X[:, npcols],
                                                      Y, cv=10)
        cvdata.append(scores)
        cvdata_avg.append(numpy.mean(scores))
        #print score
        cvlabels.append(str(len(cols)))
        #print ' ', cols, coef
        cols = cols[:maxi] + cols[maxi + 1:]

    # ge score if we always guess one value
    same_frac = numpy.sum(Y) * 1.0  / Y.shape[0]
    same_frac = max(same_frac, 1.0 - same_frac)
    print "same frac", same_frac

    for gtype in [0, 1]:
        plt.figure()
        if gtype == 0:
            plt.boxplot(numpy.array(cvdata).T)
            plt.xticks(numpy.arange(len(cvlabels)) + 1, 
                       cvlabels)
        else:
            plt.plot(cvdata_avg)
            plt.xticks(numpy.arange(len(cvlabels)), cvlabels)
            #plt.xlim((-0.5, len(cvlabels) + 1.5))
        plt.title('random forest feature reduction')
        plt.xlabel('number of features')
        plt.ylabel('score')
        plt.axhline(same_frac, color='magenta', linestyle='--')

    plt.show()

def estimateSuccess(analysis):
    #try to estimate success from long term users
    #uses multiple metrics for success such as ( having riched dist event rate 0 in biggest_improvement in analysis OR biggest improvement rate ever seen OR spearman rank and pearsons correlations between dist event rate and rate for device use OR rate for device use with "reactgood" as a reaction)
    #inputs are measurements in analysis

    #get inputs for machine learning
    X, Y, Y2, Y3, Y4=dataForML(analysis)
    #Y2, Y4 can be used for classifier
    
    #check if features separate well using MDS
    plotMDS(X, Y3)
    
    
    

    #prepare separation into train and test datasets for machine learning
    total=len(Y)
    train_num=int(total*.7) #number of data points used for train  set
    X=numpy.array(X)
    Y=numpy.array(Y)
    Y=(Y==0.0)*1
    Y2=numpy.array(Y2) * 1.0
    Y3=numpy.array(Y3) * 1.0
    Y4=numpy.array(Y4)
    print 'X', numpy.shape(X)


    #check if features separate well using PCA
    comment='pearson correlation'
    plotPCA(X, Y3, comment)
    comment='best improvement rate'  
    plotPCA(X, Y2, comment)

    #normalize data for machine learning
    X=sklearn.preprocessing.normalize(1.0* X, axis=0)
    #print 'sample X after renom', X[:10]
    #print numpy.shape(X)

    
    setSize=train_num

    
    comment='work flag'
    print '******************'
    randomForestRegressorStudy(X,Y4, setSize, comment)

    comment='pearson correlation'
    print '******************'
    randomForestRegressorStudy(X,Y3, setSize, comment)
    
    comment='best improvement rate'
    print '******************'
    randomForestRegressorStudy(X,Y2, setSize, comment)
    
    comment='best improvement rate'
    print '******************'
    linearRegressionStudy(X, Y2, setSize, comment)
    
    comment='pearson correlation'
    print '******************'
    linearRegressionStudy(X, Y3, setSize, comment) 
    
    comment='best improvement rate'
    print '******************'
    svmRegressorStudy(X, Y2, setSize, comment)
    
    comment='pearson correlation'
    print '******************'
    svmRegressorStudy(X, Y3, setSize, comment)

    comment='work flag'
    print '******************'
    logRegressionStudy(X, Y4, setSize, comment)

    comment='best rate ever achieved'
    print '******************'
    logRegressionStudy(X, Y, setSize, comment)

    comment=''
    print '******************'
    randomForestClassifierStudy(X, Y4, setSize, comment)


    print
    print 'Y4 is true', numpy.sum(Y4), '/', Y4.shape[0]

    

    # n_features = X_train.shape[1]
    # for fnum in range(31, 2**n_features):
    #     features = [i for i in range(n_features)
    #                 if (fnum >> i) & 1 == 1]
    #     print 'run %d: %s' %(fnum, features)
    #     # Create linear regression object
    #     regr = linear_model.LinearRegression()
    #     # Train the model using the training sets
    #     regr.fit(X_train[:, features], Y2_train)
    #     Y2_pred_train = regr.predict(X_train[:, features])
    #     Y2_pred_test = regr.predict(X_test[:, features])

    #     # The coefficients
    #     print('Coefficients: %s' % regr.coef_)
    #     # The mean square error
    #     print("Residual sum of squares: %.2f"
    #       % numpy.mean((Y2_pred_test - Y2_test) ** 2))
    #     # Explained variance score: 1 is perfect prediction
    #     print('Variance score: %.2f' % regr.score(
    #         X_test[:, features], Y2_test))

    # if 1:
    #     plt.subplot(1, 2, 1)
    #     plt.plot(Y2_train, Y2_pred_train, '*')
    #     plt.subplot(1, 2, 2)
    #     plt.plot(Y2_test, Y2_pred_test, '*')
    #     plt.show()


#**************************
#clustering
def clusterInit(analysis):
    #cluster based on initial conditions
    X=[]
    for i in children.keys():
        init=children[i].initial
        if analysis[i].start_date:
            X.append([init.month_count, init.night_count, len(analysis[i].use)])

    X=numpy.array(X)
    k_means = cluster.KMeans(n_clusters=3)
    k_means.fit(X) 

    #print(k_means.labels_[::])
    y=k_means.predict(X)
    total_tr=numpy.array(X)[:,0]*numpy.array(X)[:,1]
    length= numpy.array(X)[:,2]
    plt.figure()
    plt.scatter(total_tr,length, c=y, edgecolor='face')
    plt.xlabel('initial condition')
    plt.ylabel('length of use')
    plt.title('initial condition as predictor for length of use')
    plt.colorbar()
    plt.show()

    
    X2=numpy.array(zip(total_tr, length))
    k_means = cluster.KMeans(n_clusters=3)
    k_means.fit(X2) 

    #print(k_means.labels_[::])
    y2=k_means.predict(X2)
    plt.figure()
    plt.scatter(total_tr, length, c=y2, edgecolor='face')
    plt.xlabel('initial condition')
    plt.ylabel('length of use')
    plt.title('initial condition combo as predictor for length of use')
    plt.colorbar()
    plt.show()

def clusterLongTermUsers(analysis):
    #try clustering longterm users
    X=[]
    Y=[]
    Y2=[]
    Y3=[]
    for cid in analysis.keys():
        ch=analysis[cid]
        if len(ch.biggest_improvement)>0:
            vector=[ch.std_time_for_bed, ch.std_event1_time, ch.std_time_to_event1, ch.avg_time_to_event1, ch.std_offset_for_device_start, ch.avg_offset_for_device_start, ch.std_device_start_time, ch.std_device_interval, ch.use_rate, len(ch.use)]

            valid=True
            for el in vector:
                if numpy.isnan(el):
                    #print 'child with id=', cid, 'invalid data point is', vector
                    valid=False

            if valid:
                if ch.corr_dev_reactgood is None or numpy.isnan(ch.corr_dev_reactgood):
                    continue
                else:
                    Y3.append(ch.corr_dev_reactgood)
                    X.append(vector)
                    Y.append(ch.biggest_improvement[2])
                    Y2.append(ch.biggest_improvement[0])
    X=numpy.array(X)
    k_means = cluster.KMeans(n_clusters=3)
    k_means.fit(X) 

    # print "label, rate change, final rate, correlation"
    # y=k_means.predict(X)
    # for i in range(len(y)):
    #     print y[i], Y[i], Y2[i], Y3[i]
               
    X2=numpy.array(zip(Y3,numpy.array(X)[:,5], numpy.array(X)[:,5]))
    k_means = cluster.KMeans(n_clusters=5)
    k_means.fit(X2) 

   
    y=k_means.predict(X2)
    plt.figure()
    plt.scatter(X2[:,0],X2[:,1], c=y, edgecolor='face')
    plt.xlabel('correlation')
    plt.ylabel('std offset')
    plt.title('long term user analysis')
    plt.colorbar()

    plt.figure()
    plt.scatter(X2[:,0],X2[:,2], c=y, edgecolor='face')
    plt.xlabel('correlation')
    plt.ylabel('mean offset')
    plt.title('long term user analysis')
    plt.colorbar()
    
    plt.show()


#######testing

children, analysis, start, end=analysis.processData()
#estLenUseFromInit(children, analysis, start)
estimateSuccess(analysis)
