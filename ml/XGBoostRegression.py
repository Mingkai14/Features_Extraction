import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import joblib
from scripts.Error import error_obj
import os


def XGBoostRegression_Process(csv_path,save_path):
    if not os.path.exists(csv_path):
        error_obj.Something_Wrong(XGBoostRegression_Process.__name__)
        return False

    if str(csv_path).split('.')[-1]!='csv':
        error_obj.Something_Wrong(XGBoostRegression_Process.__name__)
        return False

    df = pd.read_csv(csv_path)
    df.drop('ID', axis=1, inplace=True)

    table_names = list(df.columns)
    table_names = [table_name for table_name in table_names if 'aaindex' not in table_name]

    df = df.loc[:, table_names]

    Y_reg = list(df['Experimental_DDG'])
    Y_cls = list(df['Experimental_DDG_Classification'])

    df.drop(['Experimental_DDG', 'Experimental_DDG_Classification'], axis=1, inplace=True)
    fea = df.values

    sc = StandardScaler()
    fea = sc.fit_transform(fea)
    joblib.dump(sc,save_path+'XGB_sc.m')


    X = fea
    y = Y_reg


    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


    model = XGBRegressor()
    model.fit(X_train, y_train)


    y_pred = model.predict(X_test)

    #MSE
    mse = mean_squared_error(y_test, y_pred)

    #RMSE
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    #MAE
    mae = mean_absolute_error(y_test, y_pred)

    #R^2 Score
    r2 = r2_score(y_test, y_pred)

    # Pearson
    y_test = np.array(y_test).reshape((-1, 1))
    y_pred = y_pred.reshape((-1, 1))
    yy = np.concatenate([y_test, y_pred], -1)
    yy = yy.T
    corr_matrix = np.corrcoef(yy)


    print("MSE:", mse)
    print("RMSE:", rmse)
    print("MAE:", mae)
    print("R^2 Score:", r2)
    print("Pearson:", corr_matrix[0][1])

    joblib.dump(model, save_path+'XGB.m')

    return True

    # MSE: 1.6914232880192852
    # RMSE: 1.300547303260933
    # MAE: 0.9120905282581997
    # R^2 Score: 0.6851894523969102
    # Pearson: 0.8285551614502231

def XGBoostRegression_Predict(csv_path,model_path,model_name='XGB'):
    df = pd.read_csv(csv_path)
    df.drop('ID', axis=1, inplace=True)
    table_names = list(df.columns)
    table_names = [table_name for table_name in table_names if 'aaindex' not in table_name]
    df = df.loc[:, table_names]
    df.drop(['Experimental_DDG', 'Experimental_DDG_Classification'], axis=1, inplace=True)
    fea = df.values
    sc = joblib.load(model_path+model_name+'_sc.m')
    fea = sc.transform(fea)
    model = joblib.load(model_path+model_name+'.m')
    y_pred = model.predict(fea)
    print(y_pred)
