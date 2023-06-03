
# Running
#python -m streamlit run your_script.py


"""

"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image
import subprocess
import os
import base64
import pickle
from scipy.stats import mannwhitneyu
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from chembl_webresource_client.new_client import new_client
from padelpy import from_smiles, padeldescriptor


from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split

from sklearn.metrics import mean_squared_error
from sklearn.metrics import matthews_corrcoef

from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVR
from sklearn.svm import SVC
from lazypredict.Supervised import LazyRegressor


from sklearn.model_selection import cross_val_score



from qiskit.circuit.library import ZZFeatureMap
from qiskit_machine_learning.kernels import QuantumKernel
from qiskit import transpile, BasicAer, QuantumCircuit



from lazypredict.Supervised import LazyClassifier







    
try:
    targets = get_chem_data_from(input_data)
    st.dataframe(targets[['organism', 'pref_name', 'target_components','target_type']])
except Exception as e:
    st.write("***Error due unknown target, pls type again***")

    
if targets is not None:
    ## Example for covid
    organism = targets.organism[4]
    selected_target = targets.target_chembl_id[4]



    st.write(f""" ### Study for: \n
    **Organism** : {organism} 
    **Target** : {selected_target}
    """)



    st.write(f"## Activity of {input_data} by IC50")
    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")

    df_raw = pd.DataFrame.from_dict(res)
    st.dataframe(df_raw)


    st.header("EDA")
    st.pyplot(heatmap_nulls(df_raw))


    ##DATA PREPROCESSING
    ignore = ['ligand_efficiency', 'pchembl_value', 'units'] #Ignore columns
    df_nulls = df_raw.columns[df_raw.isnull().any()] 



    df = (
        df_raw
        .isnull()
        .sum()
    )
    names = ["column_name", "null_amount"]
    df_ = pd.DataFrame({'column_names' : df_raw.columns, 'null amounts': df})

    # bar
    st.pyplot(barplot_nulls(df))

    ignore = ['ligand_efficiency', 'pchembl_value', 'units'] #Ignore columns 
    df_nulls = df_raw.columns[df_raw.isnull().any()] # find null columns.
    unwanted_columns = list(df_nulls) # Unwanted columns



    for column_to_ignore in ignore:
        unwanted_columns.remove(column_to_ignore)


    # DATA PREPROCESSING
    df_raw.drop(columns=unwanted_columns, inplace=True)

    df_transform_ = df_raw[df_raw.standard_value.notna()]
    df_transform_

    df_transform_0 = classifier(df_transform_)
    df_transform_1 = df_transform_0

    df_transform_2 = lipinksy_descriptors(df_transform_1.canonical_smiles)
    df_transform_3 =pd.concat([df_transform_0, df_transform_2], axis=1)
    df_transform_3["standard_value"] = df_transform_3["standard_value"].astype(float)
    df_transform_4 = norm_values(df_transform_3)
    df_transform_5 = pIC50(df_transform_4)
    df_transform_6 = df_transform_5[df_transform_5['bioactivity_class'] !='intermediate']


    st.header("Stats for descriptors")
    st.write("### Mannwhitney Test")
    descriptors = ['standard_value_pIC50', 'desc_molecular_weight','desc_LogP' ,'desc_hdonors','desc_hacceptors'] 
    st.dataframe(mannwhitney_test(df_transform_6,descriptors))

    st.write("### Descriptors for paDEL ")
    df_transform_6[['canonical_smiles','molecule_chembl_id']]

    padeldescriptor(mol_dir='./molecule.smi',
                d_file='PubchemFingerprinter.csv',
                #descriptortypes = 'PubchemFingerprinter.xml',
                detectaromaticity=True,
                standardizenitro=True,
                standardizetautomers=True,
                threads=2,
                removesalt=True,
                log=True,
                fingerprints=True)
    

    descriptors = pd.read_csv('./PubchemFingerprinter.csv')

    X = descriptors.drop('Name', axis=1)
    Y = df_transform_6.bioactivity_class
    #X = remove_low_variance(X, threshold=0.1)
    svm_classic = """
    
    from qiskit.circuit.library import ZZFeatureMap
    from qiskit_machine_learning.kernels import QuantumKernel
    from qiskit import transpile, BasicAer, QuantumCircuit

    # Creamos el Simulador
    backend = BasicAer.get_backend("qasm_simulator")
    shots =1024

    # Creamos el Feature Map

    dimension =2 
    feature_map =  ZZFeatureMap(dimension, reps=1)
    kernel = QuantumKernel(feature_map, quantum_instance=backend)


    # Fit an predict
    qsvm = SVC(kernel=kernel.evaluate)
    model.fit(X_train, y_train)
    model.predict(X_test)
    
    """

    manual_kernel ="""
    # Creamos una función que nos sirva para computarlo, nos retornará el valor del kernel, manualmente

    def evaluate_kernel(x_i,x_j):
        circuit = QuantumCircuit(dimension)#creo un circuito a partir de las dimensiones
        circuit.compose(feature_map.assign_parameters(x_i), inplace=True)#añado al circuito: al objeto feature map usando el metodo de asignar parámetros le doy el vector y su inversa
        circuit.compose(feature_map.assign_parameters(x_j).inverse(), inplace=True)
        circuit.measure_all()#medimos todo el circuito

        transpiled = transpile(circuit, backend)# Transpilamos el circuito
        counts =  backend.run(transpiled, shots=shots).result().get_counts() # Tomamos mediciones cantidad de veces que tenemos el estado 0

        return counts.get("0"*dimension, 0)/shots

    """

    st.header(" Quantum Support Vector Machine")
    st.code(body = svm_classic, language="python")



    # spliting data
    #X_train,X_test, y_train, y_test = train_test_split(X,Y, test_size=0.2)


    # Calsical models
    #model = SVC()
    #model.fit(X_train, y_train)