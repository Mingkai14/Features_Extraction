from ml.XGBoostRegression import XGBoostRegression_Process
from ml.SVR import SVR_Process
from ml.RFRegression import RFRegression_Process
from ml.LinearRegression import LinearRegression_Process
from ml.DTRegression import DTRegression_Process
from ml.AdaBoostRegressor import AdaBoostRegressor_Process
from ml.MLPRegression import MLPRegression_Process

from scripts.Global_Value import Features_Table_Path,Features_Table_Name

import argparse
from scripts.Error import error_obj
import os


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Input arguments')

    parser.add_argument('--model_type', type=str, default='XGB')
    parser.add_argument('--model_saving_path', type=str, default='./models/')

    args = parser.parse_args()
    if args.model_type not in ['XGB','SVR','RF','LR','DTR','ABR','MLP']:
        error_obj.Something_Wrong(__name__,'There are some wrongs in arguments')
        exit(1)

    if not os.path.exists(args.model_saving_path):
        error_obj.Something_Wrong(__name__, 'Saving path is not existed')
        exit(1)

    if args.model_type=='XGB':
        if not XGBoostRegression_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'XGB trainning failed, feature dataset may be wrong')
            exit(1)
    elif args.model_type=='SVR':
        if not SVR_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'SVR trainning failed, feature dataset may be wrong')
            exit(1)
    elif args.model_type=='RF':
        if not RFRegression_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'RF trainning failed, feature dataset may be wrong')
            exit(1)
    elif args.model_type=='LR':
        if not LinearRegression_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'LR trainning failed, feature dataset may be wrong')
            exit(1)
    elif args.model_type=='DTR':
        if not DTRegression_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'DTR trainning failed, feature dataset may be wrong')
            exit(1)
    elif args.model_type=='ABR':
        if not AdaBoostRegressor_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'ABR trainning failed, feature dataset may be wrong')
            exit(1)
    elif args.model_type=='MLP':
        if not MLPRegression_Process(Features_Table_Path + Features_Table_Name, args.model_saving_path):
            error_obj.Something_Wrong(__name__,'MLP trainning failed, feature dataset may be wrong')
            exit(1)


