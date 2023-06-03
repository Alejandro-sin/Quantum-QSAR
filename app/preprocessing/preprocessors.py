'''
Module for Preprocessing datasets for Pubchem


'''

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from sklearn.feature_selection import VarianceThreshold



## Molecular Fingerprints, Data enrichment Utilities
def classifier(data):
    '''

    Recive dataframe with a column call "standard_value" for operations.

    '''
    bioactivity_compunds =[]
    for standard_value in data.standard_value:
        standard_value = float(standard_value)
        if standard_value <=1000:
            bioactivity_compunds.append("activate")
        elif standard_value >=10000:
            bioactivity_compunds.append("inactivate")
        else:
            bioactivity_compunds.append("intermediate")
    data['bioactivity_class'] = bioactivity_compunds
    return data




## Function for calculation Lipinksy Descriptors.
def lipinksy_descriptors(smiles, verbose = False):
    i = 0
    descriptors = []
    moldata =[Chem.MolFromSmiles(molecule) for molecule in smiles]
    for mol in moldata:
        desc_molecular_weight = Descriptors.MolWt(mol)
        desc_logp = Descriptors.MolLogP(mol)
        desc_hdonors = Descriptors.NumHDonors(mol)
        desc_hacceptors = Descriptors.NumHAcceptors(mol)

        descriptors.append({
            'smile': smiles[i],
            'desc_molecular_weight': desc_molecular_weight,
            'desc_LogP': desc_logp,
            'desc_hdonors': desc_hdonors,
            'desc_hacceptors': desc_hacceptors 
        })
        i +=1

    return pd.DataFrame(data=descriptors)



def pIC50(data):
    '''
    Construct the pIC50 column in the dataframe.
    '''
    pIC50 = []
    for datapoint in data['standard_value_norm']:
        molar = datapoint * (10**-9)
        molar = -np.log10(molar)
        pIC50.append(molar)

    data['standard_value_pIC50'] =  pIC50
    return data



## Normalization function

def norm_values(data):
    '''
    Normalizes the logaritm value for returning only positive values.
    '''
    norm = []
    for datapoint in data["standard_value"]:
        if datapoint > 10**9:
            datapoint = 10**8
        norm.append(datapoint)

    data["standard_value_norm"] = norm
   #x = data.drop("standard_value", 1)
    return data




## Remove the low variance
def remove_low_variance(input_data, threshold=0.1):
    selection = VarianceThreshold(threshold)
    selection.fit(input_data)
    return input_data[input_data.columns[selection.get_support(indices=True)]]




