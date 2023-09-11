import dill
import pandas as pd
import joblib


def load_pkl(filepath):  # load model of pkl format from filepath, and return data (model)
    with open(filepath, "rb") as fr:
        data = dill.load(fr, encoding="utf-8")
    print(f"[{filepath}] data loading...")
    return data


def XGBoostRegression_Predict(csv_path,model_path,output_path):
    df = pd.read_csv(csv_path)
    ids=list(df['ID'])
    df.drop('ID', axis=1, inplace=True)

    rfe_infos = pd.read_excel(model_path+'rfe_infos.xlsx')
    X_cols = rfe_infos[rfe_infos["ranking"] == 1]["feature_names"].tolist()
    table_names = list(df.columns)
    table_names = [table_name for table_name in table_names if table_name in X_cols]
    df = df.loc[:, table_names]
    # df.drop(['Experimental_DDG', 'Experimental_DDG_Classification'], axis=1, inplace=True)
    fea = df.values
    sc = load_pkl(model_path+'sc.pkl')
    fea = sc.transform(fea)
    model = load_pkl(model_path+'predictor.pkl')
    y_pred = model.predict(fea)
    res_dict={}
    for i in range(len(y_pred)):
        res_dict[ids[i]]=y_pred[i]
    sorted_list=sorted(res_dict.items(), key=lambda x: x[1])
    print(sorted_list)
    with open(f'{output_path}/sorted_ddg.txt','w') as res:
        for item in sorted_list:
            res.write(f'{item[0]}:{item[1]}\n')


