
## import relevant packages ##
import pandas as pd
import sys


def read_prep_data(sample, sample_data):
    try:
        df = pd.read_table(sample_data, sep=',', header=0)
        df = df.loc[((df['Category']=='CLASSIFIER') |  (df['Category']=='COMBINED')) & (df['DataType']!='GENDER')]
        df.loc[df['DataType']=='SNV_96_PAIRWISE_SIMILARITY', 'DataType'] = 'SNV signatures'
        df.loc[df['DataType']=='GENOMIC_POSITION_SIMILARITY', 'DataType'] = 'somatic mutation pattern'
        df.loc[df['DataType']=='FEATURE', 'DataType'] = 'oncogenic drivers and \n features'
        df.loc[df['DataType']=='DNA_COMBINED', 'DataType'] =  'combined classifier'
        df.loc[df['RefCancerType']=='Uterus: Endometrium','RefCancerType'] = 'Endometrium'
        df.loc[df['RefCancerType']=='Colorectum/Appendix/SmallIntestine','RefCancerType'] = 'Lower GI tract'
        df['RefValue'] = df['RefValue']*100
        df = df[['DataType','RefCancerType','RefValue']].reset_index(drop=True)
        #
        df_spider = df.sort_values(['RefCancerType','DataType'],ascending= True).reset_index(drop=True)
        df_spider = df_spider[df_spider['DataType']=='combined classifier'][['RefCancerType','RefValue']]
        df_spider['reference'] = 80
        #
        df_sample_selection = df.loc[df['DataType']=='combined classifier']
        df_sample_selection = df_sample_selection.sort_values('RefValue',ascending= False).head(3).reset_index(drop=True)
        df_sample_selection = df_sample_selection[df_sample_selection['RefValue']>10]
        if len(df_sample_selection) == 0:
            df_sample_selection = df.loc[df['DataType']=='combined classifier']
            df_sample_selection = df_sample_selection.sort_values('RefValue',ascending= False).head(1).reset_index(drop=True)
        df_sample_selection['RefValue']= df_sample_selection['RefValue'].round(1)
        df_sample_selection = df_sample_selection.reset_index()
        df_sample_selection["index"] = df_sample_selection["index"]+1
        df_sample_selection["index"] = df_sample_selection["index"].astype(str)+" - "+df_sample_selection["RefCancerType"]
        df_sample_selection['index'] = df_sample_selection['index'].astype(str) + ' (likelihood=' + df_sample_selection['RefValue'].astype(str) + '%)'
        df_sample_selection = df_sample_selection[['index','RefCancerType']].reset_index(drop=True)
        df_bars = df_sample_selection.merge(df, on='RefCancerType', how='inner').sort_values(['index','DataType'],ascending= True).reset_index(drop=True)
        df_bars['RefCancerType']=df_bars['index']
        df_bars = df_bars.loc[df_bars['DataType']!='combined classifier'][['RefCancerType','DataType','RefValue']]
        #
        df = None
        df_sample_selection = None
        return df_spider, df_bars
    except:
        print('[ERROR] the sample data provided is not existing / is not in the correct format. No output files generated. CUPPA-chart will end.')
        sys.exit(1)

