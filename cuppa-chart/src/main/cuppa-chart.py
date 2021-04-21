'''USAGE:  cuppa_chart.py -sample {sampleId} -sample_data {path}/{sampleId}.cup.data.csv -output_dir {path}/{output_dir}'''

## import relevant packages ##

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import seaborn as sns
import sys
from datetime import datetime
import os

## set figure settings ##

SMALL_SIZE = 11
MEDIUM_SIZE = 11
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

## read inputs ##

in_arr = sys.argv
if '-sample' not in in_arr  or '-sample_data' not in in_arr  or '-output_dir' not in in_arr:
    print (__doc__)
    sys.exit('[ERROR] INPUT INCORRECT - not all inputs | inputs are wrongly provided. Please check USAGE above.')
else:
    sample = in_arr[in_arr.index('-sample') + 1]
    sample_data = in_arr[in_arr.index('-sample_data') + 1]
    output_dir = in_arr[in_arr.index('-output_dir') + 1]
    if output_dir.endswith('/'):
        output_dir = output_dir
    else:
        output_dir = output_dir + '/'

if sample not in sample_data:
    print (__doc__)
    sys.exit('[ERROR] INPUT INCORRECT - the sample name is different from the name on the sample data provided. Please check USAGE above.')



def main(sample: str, sample_data: str, output_dir: str):

    ## start CUPPA-chart ##
    print("CUPPA chart and conclusion generation for " + sample)
    print("Sample input: " + sample_data)
    print("Sample output: " + output_dir)
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            sys.exit('[ERROR] Output_dir dies not exist but can also not be made. No output files generated. CUPPA-chart will end.')

    ## prepare data ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - preparing sample data for " + sample)
    prepare_data(sample, sample_data)

    ## create base chart ##
    create_base_chart()

    ## add spider plot ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - creating the basis of the chart (spider plot)")
    add_spider_plot(df_spider)

    ## add barchart(s) ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - adding the relevant barcharts to the chart")
    add_barcharts(df_bars)

    ## add conclusion to figure ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - adding the final conclusion to the chart")
    add_conclusion(df_spider, df_bars)

    ## create & save conclusion file ##
    create_conclusion_file(sample, df_spider, df_bars, output_dir)
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - creating the 1st output (file with final conclusion): " + output_dir + sample + ".cuppa.conclusion.txt")

    ## create & save chart file ##
    create_chart_file(sample, df_spider, df_bars, output_dir)
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - creating the 2nd output (chart): " + output_dir + sample + ".cuppa.chart.png")

    ## CUPPA-chart complete ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [INFO] - CUPPA chart and conclusion generation for " + sample + " complete.")


def prepare_data(sample: str, sample_data: str):
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
        return(df_spider, df_bars)
    except:
        sys.exit('[ERROR] the sample data provided is not existing. No output files generated. CUPPA-chart will end.')


def create_base_chart():
    fig = plt.figure(constrained_layout=True, figsize=(13,6.8))
    gs = fig.add_gridspec(3, 44)
    return(fig, gs)


def add_spider_plot(df_spider):
    try:
        ax1 = plt.subplot(gs[:, :34], polar=True)
        ax1.set_theta_offset(pi / 2)
        ax1.set_theta_direction(-1)
        categories=list(df_spider['RefCancerType'])
        N = len(categories)
        angles = [n / float(N) * 2 * pi for n in range(N)]
        angles += angles[:1]
        plt.xticks(angles[:-1], categories)
        for label, angle in zip(ax1.get_xticklabels(), angles):
            if angle in (0, np.pi):
                label.set_horizontalalignment('center')
            elif 0 < angle < np.pi:
                label.set_horizontalalignment('left')
            else:
                label.set_horizontalalignment('right')
        ax1.set_rlabel_position(90)
        plt.yticks([0,20,40,60,80,100], ["0","20",'40','60','80','100'], color="#505050", size=10)
        plt.ylim(-6,100)

        # Add confidence
        thetas = np.linspace(0,2*np.pi,100)
        ax1.fill(thetas, [100 for i in thetas], color = "#d9f2e6", alpha = 0.4)
        ax1.fill(thetas, [80 for i in thetas], color = "#FFFFFF", alpha = 1)

        # Add plot
        values=df_spider['RefValue'].values.flatten().tolist()
        values += values[:1]
        ax1.plot(angles, values, linewidth=3, linestyle='solid', color='#0059b3', alpha=0.4, marker='.', markersize=8)
        for i in range(len(values)):
            if values[i] >= 80:
                ax1.plot(angles[i], values[i], linewidth=3, linestyle='solid', color='#00008B', marker='.', markersize=8)
                ax1.plot(angles[i-1:i+2], values[i-1:i+2], linewidth=3, linestyle='solid', color='#00008B')
        ax1.fill(angles, values, color='#0059b3', alpha=0.3)
    except:
        sys.exit('[ERROR] the basis of the chart could not be made. No output files generated. CUPPA-chart will end.')



def add_barcharts(df_bars):
    try:
        # barchart 1
        ax2 = fig.add_subplot(gs[0, 34:44])
        sns.barplot(x='RefValue',y='DataType', data=df_bars[df_bars['RefCancerType']==sorted(df_bars['RefCancerType'].unique())[0]], color="#0059b3", alpha=0.4, ax = ax2)
        ax2.set_ylabel('')
        ax2.tick_params(labelsize=9)
        ax2.set_xlabel('Per classifier likelihood (%)', style='italic')
        ax2.set_xlim(0,100)
        ax2.set_title(sorted(df_bars['RefCancerType'].unique())[0])
        labels = [item.get_text() for item in ax2.get_yticklabels()]
        locs = ax2.get_yticks()
        for loc, label in zip(locs,labels):
            if loc == 1:
                loc = loc+0.21
            ax2.text(2, loc, label, fontsize = 9)
        ax2.set_yticklabels('')
        plt.text(85, -1.25, 'Likelihood based on:', horizontalalignment='right', fontsize=14)


        if ((len(df_bars['RefCancerType'].unique())) == 1) | (df_spider['RefValue'].max() >= 80):
            ax3 = fig.add_subplot(gs[1, 34:44])
            ax3.axes.get_yaxis().set_visible(False)
            ax3.axes.get_xaxis().set_visible(False)
            ax3.set_frame_on(False)
            ax4 = fig.add_subplot(gs[2, 34:44])
            ax4.set_frame_on(False)
            ax4.axes.get_yaxis().set_visible(False)
            ax4.axes.get_xaxis().set_visible(False)
        else:
            # if applicable, add barchart 2
            ax2.set_xlabel('')
            ax3 = fig.add_subplot(gs[1, 34:44])
            sns.barplot(x='RefValue',y='DataType', data=df_bars[df_bars['RefCancerType']==sorted(df_bars['RefCancerType'].unique())[1]], color="#0059b3", alpha=0.4, ax = ax3)
            ax3.set_ylabel('')
            ax3.tick_params(labelsize=9)
            ax3.set_xlabel('Per classifier likelihood (%)', style='italic')
            ax3.set_xlim(0,100)
            ax3.set_title(sorted(df_bars['RefCancerType'].unique())[1])
            labels = [item.get_text() for item in ax3.get_yticklabels()]
            locs = ax3.get_yticks()
            for loc, label in zip(locs,labels):
                if loc == 1:
                    loc = loc+0.21
                ax3.text(2, loc, label, fontsize = 9)
            ax3.set_yticklabels('')
            if (len(df_bars['RefCancerType'].unique())) == 2:
                ax4 = fig.add_subplot(gs[2, 34:44])
                ax4.set_frame_on(False)
                ax4.axes.get_yaxis().set_visible(False)
                ax4.axes.get_xaxis().set_visible(False)
            else:
                # if applicable, add barchart 3
                ax3.set_xlabel('')
                ax4 = fig.add_subplot(gs[2, 34:44])
                sns.barplot(x='RefValue',y='DataType', data=df_bars[df_bars['RefCancerType']==sorted(df_bars['RefCancerType'].unique())[2]], color="#0059b3", alpha=0.4, ax = ax4)
                ax4.set_ylabel('')
                ax4.tick_params(labelsize=9)
                ax4.set_xlabel('Per classifier likelihood (%)', style='italic')
                ax4.set_xlim(0,100)
                ax4.set_title(sorted(df_bars['RefCancerType'].unique())[2])
                labels = [item.get_text() for item in ax4.get_yticklabels()]
                locs = ax4.get_yticks()
                for loc, label in zip(locs,labels):
                    if loc == 1:
                        loc = loc+0.21
                    ax4.text(2, loc, label, fontsize = 9)
                ax4.set_yticklabels('')
    except:
        sys.exit('[ERROR] the barcharts could not be added to the chart. No output files generated. CUPPA-chart will end.')



def add_conclusion(df_spider, df_bars):
    try:

        if (df_spider['RefValue'].max() >= 80):
            fig.suptitle('Molecular tissue of origin - ' + sorted(df_bars['RefCancerType'].unique())[0].lstrip('1 -'), fontsize=16, horizontalalignment='center')
        else:
            fig.suptitle('Molecular tissue of origin - results inconclusive', fontsize=16, horizontalalignment='center')
    except:
        sys.exit('[ERROR] the final conclusion could not be added to the chart. No output files generated. CUPPA-chart will end.')


def create_conclusion_file(sample, df_spider, df_bars, output_dir):
    try:
        conclusion_file = open(output_dir + sample + '.cuppa.conclusion.txt',"w")
        if (df_spider['RefValue'].max() >= 80):
            conclusion_file.write('Molecular tissue of origin - ' + sorted(df_bars['RefCancerType'].unique())[0].lstrip('1 -'))
        else:
            conclusion_file.write('Molecular tissue of origin - results inconclusive')
        conclusion_file.close()
    except:
        sys.exit('[ERROR] the final conclusion could not be written to a file. No output files generated. CUPPA-chart will end.')


def create_chart_file(sample, df_spider, df_bars, output_dir):
    try:
        fig.savefig(output_dir + sample + '.cuppa.chart.png', format='png', dpi=1200, bbox_inches='tight')
    except:
        os.remove(output_dir + sample + '.cuppa.conclusion.txt')
        sys.exit('[ERROR] the final chart could not be saved as a file. No output files generated. CUPPA-chart will end.')

