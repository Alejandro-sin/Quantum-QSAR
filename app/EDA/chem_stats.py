'''
Module for make and Chemical EDA of the datasets


'''
import pandas as pd
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt




def mannwhitney_test(dataframe, descriptor=None, file_csv=False):
    _doc = '''

    Dataframe: Its the target data containing the bioactivity class
    Descriptor: The "feature" or atributte that we want to compare, in our case is the descriptor
        such "pIC50" or "Molecular Weight"

    '''
    if isinstance(descriptor, str) and not None:
        active = dataframe[dataframe['bioactivity_class'] =='activate']
        inactive = dataframe[dataframe['bioactivity_class']=='inactivate']
        stats, p = mannwhitneyu(active[descriptor], inactive[descriptor])
        # Interpretation
        alpha = 0.5
        if p > alpha:
            interpretation = "Same distribution. Fail to Reject H0"
        else:
            interpretation = "Different distribution, reject H0"

        results = pd.DataFrame({
            'Descriptor': descriptor,
            'Statistics': stats,
            'p-value': p,
            'alpha': alpha,
            'Interpretation' : interpretation
        }, index =[0])
        if file_csv:
            filename =  "mannwhitneyu_results_" + descriptor + ".csv"
            path = None
            results.to_csv(filename)
        return results
    

    elif isinstance(descriptor, list) and len(descriptor) > 1:
        active = dataframe[dataframe['bioactivity_class'] =='activate']
        inactive = dataframe[dataframe['bioactivity_class']=='inactivate']

        stack_statistics = pd.DataFrame(columns=['Descriptor','Statistics','p-value','alpha','Interpretation'])
        
        for element in descriptor:
            stats, p = mannwhitneyu(active[element], inactive[element])

            alpha = 0.5
            if p > alpha:
                interpretation = "Same distribution. Fail to Reject H0"
            else:
                interpretation = "Different distribution, reject H0"

            row = pd.DataFrame({
                'Descriptor': element,
                'Statistics': stats,
                'p-value': p,
                'alpha': [alpha],
                'Interpretation': interpretation
             })

            stack_statistics = pd.concat([stack_statistics, row], ignore_index=True)

            
        return stack_statistics
    else:
        return print(_doc)





# Visualizations for
def heatmap_nulls(df_raw):
    (
    df_raw
    .isnull()
    .transpose()
    .pipe(
        lambda df: (
            sns.heatmap(data=df)
        )
    )
    )



##
def barplot_nulls(df):
    plt.figure(figsize=(15, 12))
    ax = sns.barplot(x=df.index, y=df.values,width=0.2)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    
