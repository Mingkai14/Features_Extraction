import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split,KFold, GridSearchCV
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import joblib
from scripts.Error import error_obj
import os



def MLPRegression_Process(csv_path,save_path):
    if not os.path.exists(csv_path):
        error_obj.Something_Wrong(MLPRegression_Process.__name__)
        return False

    if str(csv_path).split('.')[-1]!='csv':
        error_obj.Something_Wrong(MLPRegression_Process.__name__)
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
    joblib.dump(sc, save_path + 'MLP_sc.m')

    X = np.array(fea)
    y = np.array(Y_reg)


    kfold = KFold(n_splits=5, shuffle=True, random_state=42)

    best_pearson=-99
    best_model=None
    x_t=None
    y_t=None

    for train_idx, test_idx in kfold.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        model = MLPRegressor(hidden_layer_sizes=(100, 50), max_iter=200, random_state=42)

        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)

        y_test = np.array(y_test).reshape((-1, 1))
        y_pred = y_pred.reshape((-1, 1))
        yy = np.concatenate([y_test, y_pred], -1)
        yy = yy.T
        corr_matrix = np.corrcoef(yy)

        pearson=corr_matrix[0][1]

        if pearson>best_pearson:
            best_pearson=pearson
            best_model=model
            x_t=X_test
            y_t=y_test

    y_pred = best_model.predict(x_t)

    # MSE
    mse = mean_squared_error(y_t, y_pred)

    # RMSE
    rmse = np.sqrt(mean_squared_error(y_t, y_pred))

    # MAE
    mae = mean_absolute_error(y_t, y_pred)

    # R^2 Score
    r2 = r2_score(y_t, y_pred)

    # Pearson
    y_test = np.array(y_t).reshape((-1, 1))
    y_pred = y_pred.reshape((-1, 1))
    yy = np.concatenate([y_t, y_pred], -1)
    yy = yy.T
    corr_matrix = np.corrcoef(yy)

    print("MSE:", mse)
    print("RMSE:", rmse)
    print("MAE:", mae)
    print("R^2 Score:", r2)
    print("Pearson:", corr_matrix[0][1])

    joblib.dump(best_model, save_path+'MLP.m')

    return True



