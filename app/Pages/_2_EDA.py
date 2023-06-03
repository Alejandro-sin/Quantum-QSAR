'''

Notas, Streamlit tiene su propia forma de escribir los números par UQE LAS COSAS A LA DERECHA EN EL sidebar queden de manera consecuente.


Preguntar, como hago cuando 


'''

from eda import chem_stats
import streamlit as st

from pages._1_Target_Search import df_raw



st.pyplot(chem_stats.heatmap_nulls(df_raw))

