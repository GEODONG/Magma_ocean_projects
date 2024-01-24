# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:45:49 2022

"""

import numpy as np
import sklearn.metrics as sm
import matplotlib.pyplot as plt
import random

from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesRegressor

import sklearn.pipeline as pl
import sklearn.linear_model as lm
import sklearn.preprocessing as sp
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures


#Train the Model
def model_train(Data_Train_X,Data_Train_Y,Data_Test_X,Data_Test_Y,Regress_Mode):
    #Polynomial_Regression
    if Regress_Mode == 0:
        polynomial_features= PolynomialFeatures(degree=2)
        Data_Train_X = polynomial_features.fit_transform(Data_Train_X)

        polynomial_features= PolynomialFeatures(degree=2)
        Data_Test_X = polynomial_features.fit_transform(Data_Test_X)

        model = LinearRegression()

    #Extra Tree
    elif Regress_Mode == 1 :
        model = ExtraTreesRegressor(bootstrap=False, ccp_alpha=0.0, criterion='squared_error',
                            max_depth=20, max_features='auto', max_leaf_nodes=None,
                            max_samples=None, min_impurity_decrease=0.0,
                            min_impurity_split=None, min_samples_leaf=1,
                            min_samples_split=2, min_weight_fraction_leaf=0.0,
                            n_estimators=15, n_jobs=None, oob_score=False,
                            random_state=42, verbose=0, warm_start=False)
    elif Regress_Mode == 2 :
    #Random Forest
        model = RandomForestRegressor(bootstrap=True, ccp_alpha=0.0, criterion='squared_error',
                                max_depth=50, max_features=None, max_leaf_nodes=80,
                                max_samples=None, min_impurity_decrease=0.0,
                                min_samples_leaf=1,
                                min_samples_split=2, min_weight_fraction_leaf=0.0,
                                n_estimators=300, n_jobs=-1, oob_score=False,
                                random_state=20, verbose=0, warm_start=False) #min_impurity_split=None,

    model.fit(Data_Train_X,Data_Train_Y)
    pred_y = model.predict(Data_Train_X) #Train_Predict
    pred_yy = model.predict(Data_Test_X) #Test_Predict


    # impute the error parameter
    #R2
    r2_model = sm.r2_score(Data_Train_Y, pred_y) # impute TrainData's R2
    r2_test = sm.r2_score(Data_Test_Y, pred_yy)  #R2 impute TestData's R2

    print(r2_model)

    
#RMSE
    rmse_model = np.sqrt(sm.mean_squared_error(Data_Train_Y, pred_y)) #RMSE impute TrainData's RMESE
    rmse_test = np.sqrt(sm.mean_squared_error(Data_Test_Y, pred_yy)) # impute TestData's RMSE

     # Keep three decimals
    r2_model = round(r2_model,3)
    r2_test = round(r2_test,3)
    rmse_model = round(rmse_model,3)
    rmse_test = round(rmse_test,3)

    #tempr = Error Parameter Set
    temp_r =[]
    temp_r.append(r2_model)
    temp_r.append(r2_test)
    temp_r.append(rmse_model)
    temp_r.append(rmse_test)


    # plot
    binary_plot(y_train = Data_Train_Y,
                y_train_label = pred_y,
                y_test = Data_Test_Y,
                y_test_label = pred_yy,
                train_rmse = rmse_model,
                test_rmse = rmse_test,
                train_r2 = r2_model,
                test_r2 = r2_test)
    save_fig("Result_plot")



# make plot by sany He
def binary_plot(y_train,  y_train_label, y_test, y_test_label,
                train_rmse, test_rmse, train_r2, test_r2,
                text_position=[0.5, -0.075]):
    """plot the binary diagram

    :param y_train: the label of the training data set
    :param y_train_label: the prediction of the training the data set
    :param y_test: the label of the testing data set
    :param y_test_label: the prediction of the testing data set
    :param train_rmse: the RMSE score of the training data set
    :param test_rmse: the RMSE score of the testing data set
    :param train_r2: the R2 score of the training data set
    :param test_r2: the R2 score of the testing data set
    :param test_position: the coordinates of R2 text for
    """
    
    plt.figure(figsize=(6,6))
    plt.scatter(y_train, y_train_label, marker="s",c='none',edgecolors='black',
                label="RMSE(train) = {}".format(train_rmse)+'\n'+r'     $R^2$(train) = {}'.format(train_r2))
    plt.scatter(y_test, y_test_label, marker="o",
                label=" RMSE(test) = {}".format(test_rmse)+'\n'+r'      $R^2$(test) = {}'.format(test_r2))
    plt.legend(loc="lower right", fontsize=14,fancybox=False,edgecolor= 'black')
    plt.xlabel("Reference value", fontsize=20,font='Arial')#: \n Using New "+r'$\Delta G$'+' Result', fontsize=16)
    plt.ylabel("Predicted value", fontsize=20,font='Arial')
    a=[0,3]; b=[0,3]
    plt.plot(a, b,c='black')
    plt.text(x= 3, y = 3,s='1:1 line',rotation = 45,fontsize= 16,font='Arial')
    a=[3.7,4]; b=[3.7,4]
    plt.plot(a, b,c='black')
    #plt.text(text_position[0], text_position[1]+2.75,
             #r'$R^2$(train)={}'.format(train_r2),
             #fontdict={'family':'Arial','size': 16, 'color': '#000000'})
    #plt.text(text_position[0], text_position[1]+2.5,
            # r'$R^2$(test)={}'.format(test_r2),
            # fontdict={'family':'Arial','size': 16, 'color': '#000000'})
    plt.xticks(np.arange(-1, 4.1), fontsize=20)
    plt.xticks(np.arange(-1,4.1,0.1),fontsize=20,minor=True)
    plt.yticks(np.arange(-1, 4.1), fontsize=20)
    plt.yticks(np.arange(-1,4.1,0.1),fontsize=20,minor= True)
    plt.xlim(-1,4)
    plt.ylim((-1, 4))


#save figure by sany He
def save_fig(fig_id, tight_layout=True):
    '''
    Run to save automatic pictures
    
    :param fig_id: image name
    '''
    path = "./Image/"+fig_id+".pdf"
    print("Saving figure", fig_id)
    if tight_layout:
        plt.tight_layout()
    plt.savefig(path, format='pdf', dpi=300)





#input data：data_ is source_data，data_0 is "label = 0" data，data_1 is "label = 1"data,data_2 is"label = 2"data
import pandas as pd
data_ = pd.read_csv("./run_new_2023_5_13b.csv", low_memory=False)
data_0 = np.loadtxt('./_data0.txt')
data_1 = np.loadtxt('./_data1.txt')
data_2 = np.loadtxt('./_data2.txt')


#Disarrange the sample order of the data set
# np.random.shuffle(data_)
np.random.shuffle(data_0)
np.random.shuffle(data_1)
np.random.shuffle(data_2)
 # load dataset

# np.random.shuffle(dataset)


#DataMode：Stratified sampling was conducted directly (7:3)
x1 = data_0[0:int(data_0.shape[0]*0.7)] 
x2 = data_1[0:int(data_1.shape[0]*0.7)] 
x3 = data_2[0:int(data_2.shape[0]*0.7)]
Data_Train = np.vstack((x1,x2,x3))

x1 = data_0[int(data_0.shape[0]*0.7):int(data_0.shape[0])] 
x2 = data_1[int(data_1.shape[0]*0.7):int(data_1.shape[0])] 
x3 = data_2[int(data_2.shape[0]*0.7):int(data_2.shape[0])]
Data_Test = np.vstack((x1,x2,x3))

Data_Train_X = Data_Train[:,0:11]
Data_Train_Y = Data_Train[:,-1]
Data_Test_X = Data_Test[:,0:11]
Data_Test_Y = Data_Test[:,-1]

Regress_Mode = 2  #0-linear 1-extratree 2-Random forest
list_R2_RMSE = model_train(Data_Train_X, Data_Train_Y, Data_Test_X, Data_Test_Y,Regress_Mode)


# import shap
# import sklearn
# regressor = ensemble.RandomForestRegressor()
# regressor.fit(X_train, y_train);
# # Create object that can calculate shap values
# explainer = shap.TreeExplainer(model)


