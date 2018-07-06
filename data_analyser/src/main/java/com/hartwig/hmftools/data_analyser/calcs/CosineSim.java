package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Double.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CosineSim {

    private static final Logger LOGGER = LogManager.getLogger(CosineSim.class);

    private String mOutputDir;

    public static int CSSR_I1 = 0;
    public static int CSSR_I2 = 1;
    public static int CSSR_VAL = 2;

    private static double MIN_VALUE = 1e-10;

    public CosineSim(final String outputDir)
    {
        mOutputDir = outputDir;
    }

    public static double calcCSS(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length || set1.length == 0)
            return 0;

        double aaTotal = 0;
        double bbTotal = 0;
        double abTotal = 0;

        for(int i = 0; i < set1.length; ++i)
        {
            double a = set1[i];
            double b = set2[i];

            aaTotal += a*a;
            bbTotal += b*b;
            abTotal += a*b;
        }

        if(aaTotal <= 0 || bbTotal <= 0)
            return 0;

        return min(abTotal / (sqrt(aaTotal) * sqrt(bbTotal)), 1.0);
    }

    public static List<double[]> getTopCssPairs(final NmfMatrix matrix1, final NmfMatrix matrix2, double cssMatchCutoff, boolean applyExclusivity, boolean skipRepeats)
    {
        // use CSS to compare each pair of columns from the 2 sets
        // returns a list of results where the CSS value is above the cutoff, ordered by CSS descending (ie closest first)

        List<double[]> cssResults = Lists.newArrayList();

        if(matrix1.Rows != matrix2.Rows)
            return cssResults;

        // record each combination of vector comparisons

        for(int i = 0; i < matrix1.Cols; ++i) {

            double[] data1 = matrix1.getCol(i);

            int j = 0;

            if(skipRepeats)
                j = i + 1;

            for (; j < matrix2.Cols; ++j) {

                double[] data2 = matrix2.getCol(j);

                double css = CosineSim.calcCSS(data1, data2);

                if (css < cssMatchCutoff)
                    continue;

                int index = 0;
                for(;index < cssResults.size(); ++index)
                {
                    final double[] result = cssResults.get(index);
                    if(css > result[CSSR_VAL])
                        break;
                }

                double[] result = {i, j, css};
                cssResults.add(index, result);
            }
        }

        if(!applyExclusivity)
            return cssResults;

        // now find the top CSS for each sig in the first set
        for(int i = 0; i < matrix1.Cols; ++i) {

            int index = 0;
            boolean topFound = false;
            double[] topResult = null;

            while (index < cssResults.size()) {

                final double[] result = cssResults.get(index);

                if(topFound)
                {
                    if((result[CSSR_I1] == topResult[CSSR_I1] || result[CSSR_I2] == topResult[CSSR_I2]))
                    {
                        cssResults.remove(index);
                    }
                    else
                    {
                        ++index;
                    }

                    continue;
                }

                if (result[CSSR_I1] != i) {
                    ++index;
                    continue;
                }

                topResult = result;
                topFound = true;
                ++index;
            }
        }

        return cssResults;
    }

    public static void logSimilarites(final NmfMatrix matrix, double cssMatchCutoff, final String item)
    {
        // use CSS to compare each pair of values and log similar ones

        for (int i = 0; i < matrix.Cols; ++i) {

            double[] sig1 = matrix.getCol(i);

            for (int j = i+1; j < matrix.Cols; ++j) {

                double[] sig2 = matrix.getCol(j);

                double css = CosineSim.calcCSS(sig1, sig2);

                if(css > 1)
                {
                    LOGGER.warn("CSS above 1");
                }

                if (css < cssMatchCutoff)
                    continue;

                LOGGER.debug(String.format("close CSS data: %s1(%d) vs %s2(%d) css(%.4f)", item, i, item, j, css));
            }
        }
    }


    public void calcCosineSimilarities(final List<String> itemIds, final List<List<Double>> dataSets, double cutoff)
    {
        if(dataSets.isEmpty())
            return;

        LOGGER.debug("calc CSS for {} items", itemIds.size());

        try
        {
            BufferedWriter writer = DataUtils.getNewFile(mOutputDir, "css_values,csv");
            writer.write("Sample1,Sample2,CSS\n");

            double[][] dataArray = DataUtils.convertArray(dataSets, true);
            int itemCount = itemIds.size();

            for(int i = 0; i < itemCount; ++i)
            {
                double[] set1 = dataArray[i];

                for(int j = i+1; j < itemCount; ++j) {

                    double[] set2 = dataArray[j];

                    double css = calcCSS(set1, set2);

                    if(css >= cutoff)
                    {
                        // LOGGER.debug("items({} and {}) css({})", itemIds.get(i), itemIds.get(j), css);
                        writer.write(String.format("%s,%s,%.6f", itemIds.get(i), itemIds.get(j), css));
                        writer.newLine();
                    }
                }

                if(i > 0 && (i % 100) == 0)
                {
                    LOGGER.debug("processed {} items", i);
                }
            }

            if(writer != null)
                writer.close();

        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }



}
