import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import AdaBoostRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import joblib
from scripts.Error import error_obj
import os


def AdaBoostRegressor_Process(csv_path,save_path):
    if not os.path.exists(csv_path):
        error_obj.Something_Wrong(AdaBoostRegressor_Process.__name__)
        return False

    if str(csv_path).split('.')[-1]!='csv':
        error_obj.Something_Wrong(AdaBoostRegressor_Process.__name__)
        return False


    df = pd.read_csv(csv_path)  # 1:未违约；0：违约
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
    joblib.dump(sc, save_path + 'ABR_sc.m')

    X = fea
    y = Y_reg

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


    model = AdaBoostRegressor()
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    mse = mean_squared_error(y_test, y_pred)

    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    mae = mean_absolute_error(y_test, y_pred)

    r2 = r2_score(y_test, y_pred)

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

    joblib.dump(model, save_path + 'ABR.m')

    return True

# 均方误差（MSE）: 2.9321150076252303
# 均方根误差（RMSE）: 1.7123419657373438
# 平均绝对误差（MAE）: 1.3064772099229314
# 决定系数（R^2 Score）: 0.4542698225074858
# 皮尔逊相关系数: 0.6932105168696464

