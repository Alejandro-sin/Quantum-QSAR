

import streamlit as st
import pandas as pd
import time
import re



#data ingestion
from ingestion.chem_ingestors import get_chem_data_from, get_activity_data

#Preprocessing





## ST Options
st.set_option('deprecation.showPyplotGlobalUse', False)

# Using object notation
input_selectbox = st.sidebar.selectbox(
    "Select example organism or protein to study",
    ("Coronavirus", "Acetil", "Reumatoide...")
)


# Para cuando tendremos una opción personalizada para que el usuario digite.
####
# Search en PubChem databse
#st.header('Type your Target')
#input_data = str(st.text_input(''))
#if input_data is None:
#    st.write("Write something for PubChem data")



#with st.spinner("Loading ..."):
    #time.sleep(1)
df_response = get_chem_data_from(input_selectbox)
#Primer Input Data
st.header("Organism or Protein")
st.markdown("""

Select an organism or protein to search at Pubchem database.



-

_____________________________________
Notas para hacer en este front End:

-  Poner isntrucciones de manejo
-  Ampliar el tamaño de la tabla para que se peciban mejor.
-  Mover la columna Selection al principio del daframe.
-  Usar la estrategia de cahe para guardar los estados de las apps.


""")
df_chem_data = df_response[['organism', 'pref_name', 'target_type', 'target_chembl_id']]
df_chem_data['Selection'] = pd.Series([False] * len(df_chem_data))
df_edited = st.experimental_data_editor(df_chem_data) # Quiero usar la función que desaparece la data para solo mostrar la opción.
selection_to_study = df_edited.loc[df_edited["Selection"]]
organism_str =  selection_to_study['organism'].values[0]
chem_id = selection_to_study['target_chembl_id'].values[0]
p_name = selection_to_study['pref_name'].values[0]

st.write(organism_str)

## Texto de ayuda
st.text(f"""

Organism: "{organism_str}" 
Target Chem  ID: "{chem_id}"
Prefered Name: "{p_name}"


Aquí implementar la API de Bard par apoder decirle que interprete y retorne descripción técncia al respecto. O cualquier LLM que me permita retornar:

1. Descripción detallada de: ¿ Qué es el tan. ..


2. Podmeo poner la molecula quimica como se, y la descripción, es decir, mostrar la data más amigable.

""")
#organism
df_raw = get_activity_data(chem_id)
st.dataframe(df_raw)
    
df_raw
