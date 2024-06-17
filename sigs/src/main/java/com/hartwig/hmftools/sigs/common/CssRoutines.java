package com.hartwig.hmftools.sigs.common;

import static java.lang.Math.abs;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.stats.CosineSimilarity.calcCosineSim;
import static com.hartwig.hmftools.common.utils.Arrays.equalArray;
import static com.hartwig.hmftools.common.utils.VectorUtils.copyVector;
import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

public class CssRoutines
{
    public static final int CSSR_I1 = 0;
    public static final int CSSR_I2 = 1;
    public static final int CSSR_VAL = 2;

    public static int getLeastSimilarEntry(final double[] set1, final double[] set2)
    {
        if(set1.length != set2.length)
            return -1;

        double refCss = calcCosineSim(set1, set2, true);

        double[] copy = new double[set1.length];
        copyVector(set1, copy);
        int bestIndex = -1;
        double bestCss = refCss;

        for(int i = 0; i < set1.length; ++i)
        {
            if(set1[i] == 0 || set2[i] == 0)
                continue;

            copy[i] = 0;

            double css = calcCosineSim(copy, set2, true);

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

    public static List<double[]> getTopCssPairs(
            final Matrix matrix1, final Matrix matrix2, double cssMatchCutoff,
            boolean applyExclusivity, boolean skipRepeats)
    {
        return getTopCssPairs(matrix1, matrix2, cssMatchCutoff, applyExclusivity, skipRepeats, false, false);
    }

    public static List<double[]> getTopCssPairs(
            final Matrix matrix1, final Matrix matrix2, double cssMatchCutoff,
            boolean applyExclusivity, boolean skipRepeats, boolean skipAllZeros, boolean skipZeroEntries)
    {
        // use CSS to compare each pair of columns from the 2 sets
        // returns a list of results where the CSS value is above the cutoff, ordered by CSS descending (ie closest first)

        List<double[]> cssResults = Lists.newArrayList();

        if(matrix1.Rows != matrix2.Rows)
            return cssResults;

        double[] emptyData = new double[matrix1.Rows];

        // record each combination of vector comparisons
        for(int i = 0; i < matrix1.Cols; ++i) {

            /*
            if(i > 0 && (i % 100) == 0)
            {
                SIG_LOGGER.debug("processed {} items", i);
            }
            */

            double[] data1 = matrix1.getCol(i);

            if(skipAllZeros && equalArray(data1, emptyData))
                continue;

            int j = 0;

            if(skipRepeats)
                j = i + 1;

            for (; j < matrix2.Cols; ++j) {

                double[] data2 = matrix2.getCol(j);

                if(skipAllZeros && equalArray(data2, emptyData))
                    continue;

                double css = calcCosineSim(data1, data2, skipZeroEntries);

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

    public static void logSimilarites(final Matrix matrix, double cssMatchCutoff, final String item)
    {
        // use CSS to compare each pair of values and log similar ones
        for (int i = 0; i < matrix.Cols; ++i) {

            double[] sig1 = matrix.getCol(i);

            for (int j = i+1; j < matrix.Cols; ++j) {

                double[] sig2 = matrix.getCol(j);

                double css = calcCosineSim(sig1, sig2);

                if(css > 1)
                {
                    CommonUtils.SIG_LOGGER.warn("CSS above 1");
                }

                if (css < cssMatchCutoff)
                    continue;

                CommonUtils.SIG_LOGGER.debug(String.format("close CSS data: %s(%d) vs %s(%d) css(%.4f)", item, i, item, j, css));
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
            CommonUtils.SIG_LOGGER.debug(String.format("prob(%.8f) totals(test=%.1f diff=%.1f same=%.1f) degFreedom(%d)",
                    probability, testVal, diffTotal, sameTotal, degFreedom));
        }

        return probability;
    }

    public static List<double[]> getTopLogLikelihoodPairs(final Matrix matrix1, final Matrix matrix2, double probCutoff)
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
                CommonUtils.SIG_LOGGER.debug("processed {} items", i);
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


}
