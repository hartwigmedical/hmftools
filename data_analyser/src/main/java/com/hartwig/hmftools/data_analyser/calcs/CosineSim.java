package com.hartwig.hmftools.data_analyser.calcs;

import static java.lang.Double.doubleToLongBits;
import static java.lang.Double.max;
import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.lang.Math.sqrt;

import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.copyVector;
import static com.hartwig.hmftools.data_analyser.calcs.DataUtils.sumVector;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.data_analyser.types.NmfMatrix;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
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
        return calcCSS(set1, set2, false);
    }

    public static double calcCSS(final double[] set1, final double[] set2, boolean skipZeros)
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

            if(skipZeros && (a == 0 || b == 0))
                continue;

            aaTotal += a*a;
            bbTotal += b*b;
            abTotal += a*b;
        }

        if(aaTotal <= 0 || bbTotal <= 0)
            return 0;

        return min(abTotal / (sqrt(aaTotal) * sqrt(bbTotal)), 1.0);
    }

    public static double calcCSSRelative(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length || set1.length == 0)
            return 0;

        double aaTotal = 0;
        double bbTotal = 0;
        double abTotal = 0;

        // reduce each pair of values to be roughly in line with the median pair
        double total1 = sumVector(set1);
        double total2 = sumVector(set2);

        if(total1 == 0 || total2 == 0)
            return 0;

        double total = (total1 + total2) * 0.5;

        for(int i = 0; i < set1.length; ++i)
        {
            double a = set1[i];
            double b = set2[i];

            if(a == 0 || b == 0)
                continue;

            double ratioToTotal = ((a + b) * 0.5) / total;

            a /= ratioToTotal;
            b /= ratioToTotal;

            aaTotal += a*a;
            bbTotal += b*b;
            abTotal += a*b;
        }

        if(aaTotal <= 0 || bbTotal <= 0)
            return 0;

        return min(abTotal / (sqrt(aaTotal) * sqrt(bbTotal)), 1.0);
    }

    public static int getLeastSimilarEntry(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length)
            return -1;

        double refCss = calcCSS(set1, set2, true);

        double[] copy = new double[set1.length];
        copyVector(set1, copy);
        int bestIndex = -1;
        double bestCss = refCss;

        for(int i = 0; i < set1.length; ++i)
        {
            if(set1[i] == 0 || set2[i] == 0)
                continue;

            copy[i] = 0;

            double css = calcCSS(copy, set2, true);

            if(css > bestCss)
            {
                bestCss = css;
                bestIndex = i;
            }

            // restore
            copy[i] = set1[i];
        }

        return bestIndex;
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

            if(i > 0 && (i % 100) == 0)
            {
                LOGGER.debug("processed {} items", i);
            }

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

                LOGGER.debug(String.format("close CSS data: %s(%d) vs %s(%d) css(%.4f)", item, i, item, j, css));
            }
        }
    }

    public static double calcLogLikelihood(final double[] set1, final double[] set2, boolean logVerbose)
    {
        if(set1.length != set2.length || set1.length == 0)
            return 0;

        double diffTotal = 0;
        double sameTotal = 0;

        double aTotal = sumVector(set1);
        double bTotal = sumVector(set2);
        double ratio = aTotal/bTotal;

        for(int i = 0; i < set1.length; ++i)
        {
            int a = (int)set1[i];
            int b = (int)round(set2[i] * ratio);

            if(a == 0 && b == 0)
                continue;

            PoissonDistribution poisson = null;
            if(a == 0)
            {
                poisson = new PoissonDistribution(b);
                diffTotal += poisson.logProbability(a);
                sameTotal += poisson.logProbability(b);
            }
            else
            {
                poisson = new PoissonDistribution(a);
                diffTotal += poisson.logProbability(b);
                sameTotal += poisson.logProbability(a);
            }
        }

        int degFreedom = set1.length - 1;
        ChiSquaredDistribution chiSquDist = new ChiSquaredDistribution(degFreedom);

        double testVal = diffTotal - sameTotal;
        double chiInput = abs(2 * testVal);

        double probability = 1 - chiSquDist.cumulativeProbability(chiInput);

        if(logVerbose)
        {
            LOGGER.debug(String.format("prob(%.8f) totals(test=%.1f diff=%.1f same=%.1f) degFreedom(%d)",
                    probability, testVal, diffTotal, sameTotal, degFreedom));
        }

        return probability;
    }

    public static List<double[]> getTopLogLikelihoodPairs(final NmfMatrix matrix1, final NmfMatrix matrix2, double probCutoff)
    {
        // use CSS to compare each pair of columns from the 2 sets
        // returns a list of results where the CSS value is above the cutoff, ordered by CSS descending (ie closest first)

        List<double[]> lliResults = Lists.newArrayList();

        if(matrix1.Rows != matrix2.Rows)
            return lliResults;

        // record each combination of vector comparisons

        for(int i = 0; i < matrix1.Cols; ++i) {

            if(i > 0 && (i % 100) == 0)
            {
                LOGGER.debug("processed {} items", i);
            }

            double[] data1 = matrix1.getCol(i);

            int j = i + 1;

            for (; j < matrix2.Cols; ++j) {

                double[] data2 = matrix2.getCol(j);

                double llProb = calcLogLikelihood(data1, data2, false);

                if (llProb < probCutoff)
                    continue;

                int index = 0;
                for(;index < lliResults.size(); ++index)
                {
                    final double[] result = lliResults.get(index);
                    if(llProb > result[CSSR_VAL])
                        break;
                }

                double[] result = {i, j, llProb};
                lliResults.add(index, result);
            }
        }

       return lliResults;
    }

    public static int calcPoissonRangeGivenProb(int value, double requiredProb)
    {
        if(value <= 10)
            return 9;

        PoissonDistribution poisson = new PoissonDistribution(value);

        int maxIterations = 10;
        int iterations = 0;

        double initRange = 3.7 / sqrt(value); // works for requiredProb = 1e-4
        int testValue = (int)max(round(value * (1 - initRange)), 0);
        int testValueUpper = (int)max(round(value * (1 - initRange*0.5)), 0);
        int testValueLower = (int)max(round(value * (1 - initRange*2)), 0);

        double currentProb = poisson.cumulativeProbability(testValue);
        double probDiff = 0;

        while(iterations < maxIterations)
        {
            probDiff = abs(requiredProb - currentProb) / requiredProb;

            if(probDiff < 0.1)
                break;

            // if prob is too high, need to lower the test value
            if(currentProb > requiredProb)
            {
                if(testValue <= testValueLower + 1)
                    break;

                testValueUpper = testValue;
                testValue = (int)round((testValue + testValueLower) * 0.5);
            }
            else
            {
                if(testValue >= testValueUpper - 1)
                    break;

                testValueLower = testValue;
                testValue = (int)round((testValue + testValueUpper) * 0.5);
            }

            currentProb = poisson.cumulativeProbability(testValue);
            ++iterations;
        }

        if(iterations >= maxIterations)
        {
            LOGGER.warn(String.format("max iterations reached: value(%d) test(%d) prob(%.4f diff=%.4f)", value, testValue, currentProb, probDiff));
        }

        return value - testValue;
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
