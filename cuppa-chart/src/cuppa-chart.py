'''USAGE:  prepare_data.py -sample {sampleId} -sample_data {path}/{sampleId}.cup.data.csv -output_dir {path}/{output_dir}'''

try:
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
    import main.prepare_data
    import main.create_chart
except:
    print('[ERROR] Dependencies of CUPPA-chart not existing. Please check your environment. CUPPA-chart can not run.')
    exit()


def main(sample, sample_data, output_dir):

    ## start CUPPA-chart ##
    print("CUPPA chart and conclusion generation for " + sample)
    print("Sample input: " + sample_data)
    print("Sample output: " + output_dir)
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            sys.exit('[ERROR] Output_dir does not exist but can also not be made. No output files generated. CUPPA-chart will end.')

    ## prepare data ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - preparing sample data for " + sample)
    df_spider, df_bars = prepare_data(sample, sample_data)

    ## create base chart ##
    fig, gs = create_base_chart()

    ## add spider plot ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - creating the basis of the chart (spider plot)")
    fig, gs = add_spider_plot(df_spider, fig, gs)

    ## add barchart(s) ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - adding the relevant barcharts to the chart")
    fig, gs = add_barcharts(df_bars, df_spider, fig, gs)

    ## add conclusion to figure ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - adding the final conclusion to the chart")
    fig, gs = add_conclusion(df_spider, df_bars, fig, gs)

    ## create & save conclusion file ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - creating the 1st output (file with final conclusion): " + output_dir + sample + ".cuppa.conclusion.txt")
    create_conclusion_file(sample, df_spider, df_bars, output_dir)

    ## create & save chart file ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [DEBUG] - creating the 2nd output (chart): " + output_dir + sample + ".cuppa.chart.png")
    create_chart_file(sample, df_spider, df_bars, output_dir, fig, gs)

    ## CUPPA-chart complete ##
    print(datetime.now().strftime("%Y/%m/%d %H:%M:%S") + " - [INFO] - CUPPA chart and conclusion generation for " + sample + " complete.")


if __name__ == "__main__":
    # read inputs
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
    # start main function
    main(sample, sample_data, output_dir)
