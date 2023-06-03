'''
Module for get the data from DB Pubchem


'''

import pandas as pd
from chembl_webresource_client.new_client import new_client



def get_chem_data_from(organism,target_protein=None):
    target =  new_client.target # Creamos el cliente
    target_query = target.search(organism)
    targets = pd.DataFrame.from_dict(target_query)
    return targets


def get_activity_data(selected_target):
    activity = new_client.activity
    res = activity\
        .filter(target_chembl_id=selected_target)\
        .filter(standard_type="IC50")
    df_raw = pd.DataFrame.from_dict(res)
    return df_raw