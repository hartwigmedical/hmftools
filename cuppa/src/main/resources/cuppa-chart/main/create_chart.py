
## import relevant packages ##
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from math import pi
import seaborn as sns
import sys
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


def create_base_chart():
    fig = plt.figure(constrained_layout=True, figsize=(13,6.8))
    gs = fig.add_gridspec(3, 44)
    return fig, gs


def add_spider_plot(df_spider, fig, gs):
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
        return fig, gs
    except:
        print('[ERROR] the basis of the chart could not be made. No output files generated. CUPPA-chart will end.')
        sys.exit(1)



def add_barcharts(df_bars, df_spider, fig, gs):
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
        return fig, gs
    except:
        print('[ERROR] the barcharts could not be added to the chart. No output files generated. CUPPA-chart will end.')
        sys.exit(1)



def add_conclusion(df_spider, df_bars, fig, gs):
    try:
        if (df_spider['RefValue'].max() >= 80):
            fig.suptitle('Molecular tissue of origin - ' + sorted(df_bars['RefCancerType'].unique())[0].lstrip('1 -'), fontsize=16, horizontalalignment='center')
        else:
            fig.suptitle('Molecular tissue of origin - results inconclusive', fontsize=16, horizontalalignment='center')
        return fig, gs
    except:
        print('[ERROR] the final conclusion could not be added to the chart. No output files generated. CUPPA-chart will end.')
        sys.exit(1)


def create_conclusion_file(sample, df_spider, df_bars, output_dir):
    try:
        conclusion_file = open(output_dir + sample + '.cuppa.conclusion.txt',"w")
        if (df_spider['RefValue'].max() >= 80):
            conclusion_file.write('Molecular tissue of origin - ' + sorted(df_bars['RefCancerType'].unique())[0].lstrip('1 -'))
        else:
            conclusion_file.write('Molecular tissue of origin - results inconclusive')
        conclusion_file.close()
    except:
        print('[ERROR] the final conclusion could not be written to a file. No output files generated. CUPPA-chart will end.')
        sys.exit(1)


def create_chart_file(sample, df_spider, df_bars, output_dir, fig, gs):
    try:
        fig.savefig(output_dir + sample + '.cuppa.chart.png', format='png', dpi=1200, bbox_inches='tight')
    except:
        os.remove(output_dir + sample + '.cuppa.conclusion.txt')
        print('[ERROR] the final chart could not be saved as a file. No output files generated. CUPPA-chart will end.')
        sys.exit(1)



