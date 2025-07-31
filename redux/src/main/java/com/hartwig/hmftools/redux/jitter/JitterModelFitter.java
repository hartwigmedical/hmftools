package com.hartwig.hmftools.redux.jitter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.redux.JitterCountsTable;
import com.hartwig.hmftools.common.redux.JitterModelParams;
import com.hartwig.hmftools.common.redux.JitterTableRow;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.jetbrains.annotations.Nullable;

public class JitterModelFitter
{
    public static class ScaleSkew
    {
        public final double scale;
        public final double skew;
        public final double loss;

        private ScaleSkew(final double scale, final double skew, final double loss)
        {
            this.scale = scale;
            this.skew = skew;
            this.loss = loss;
        }

        public static ScaleSkew of(final double scale, final double skew, final double loss)
        {
            return new ScaleSkew(scale, skew, loss);
        }
    }

    public static class FittingDefaults
    {
        public final double scale_four;
        public final double scale_five;
        public final double scale_six;
        public final double slope;

        private FittingDefaults(final double scale_four, final double scale_five, final double scale_six, final double slope)
        {
            this.scale_four = scale_four;
            this.scale_five = scale_five;
            this.scale_six = scale_six;
            this.slope = slope;
        }
    }

    public FittingDefaults defaults(int repeatUnitLength)
    {
        switch(repeatUnitLength)
        {
            case 1:
                return new FittingDefaults(0.1, 0.15, 0.2, 0.06);
            case 2:
                return new FittingDefaults(0.15, 0.2, 0.25, 0.06);
            default:
                return new FittingDefaults(0.2, 0.25, 0.3, 0.06);
        }
    }

    private final JitterCountsTable mStatsTable;
    private final Map<Integer, ScaleSkew> msRepeatFittedScaleSkew = new HashMap<>();

    private JitterModelParams mJitterModelParams = null;

    public JitterModelParams getJitterModelParams()
    {
        return mJitterModelParams;
    }

    public JitterModelFitter(final JitterCountsTable statsTable)
    {
        mStatsTable = statsTable;
    }

    public void performFit()
    {
        int repeatUnitLength = mStatsTable.repeatUnitLength();
        FittingDefaults relevantDefaults = defaults(repeatUnitLength);
        msRepeatFittedScaleSkew.clear();
        mJitterModelParams = null;
        int totalShortenedCount = 0;
        int totalLengthenedCount = 0;

        for(int numRepeats = 4; numRepeats <= 15; ++numRepeats)
        {
            double defaultScale = Double.NaN;
            switch(numRepeats)
            {
                case 4:
                    defaultScale = relevantDefaults.scale_four;
                    break;
                case 5:
                    defaultScale = relevantDefaults.scale_five;
                    break;
                case 6:
                    defaultScale = relevantDefaults.scale_six;
                    break;
            }
            JitterTableRow row = mStatsTable.getRow(numRepeats);
            if(row != null)
            {
                for(Map.Entry<Integer, Integer> entry : row.jitterCounts().entrySet())
                {
                    if(entry.getKey() < 0) { totalShortenedCount += entry.getValue(); }
                    if(entry.getKey() > 0) { totalLengthenedCount += entry.getValue(); }
                }
                ScaleSkew bestScaleSkew = gridSearch(numRepeats, row, defaultScale);
                msRepeatFittedScaleSkew.put(numRepeats, bestScaleSkew);
            }
        }

        double readCountSum7_15 = 0.0;
        int numRepeatCountsOverMin = 0;
        for(int numRepeats = 7; numRepeats <= 15; ++numRepeats)
        {
            int readCount = mStatsTable.getReadCount(numRepeats);
            readCountSum7_15 += readCount;
            if(readCount >= 20000) { numRepeatCountsOverMin++; }
        }

        double scaleFitGradient = Double.NaN;
        double scaleFitIntercept = Double.NaN;
        double microsatelliteSkew = Double.NaN;

        if(numRepeatCountsOverMin > 1)
        {
            // next we do a weighted least square fitting linear regression through the repeat lengths
            // 7-15 to find the scale
            // https://stackoverflow.com/questions/5684282/weighted-linear-regression-in-java

            List<Double> xUnweighted = new ArrayList<>();
            List<Double> yUnweighted = new ArrayList<>();
            List<Double> weights = new ArrayList<>();
            List<Double> skews = new ArrayList<>();

            for(int numRepeats = 7; numRepeats <= 15; ++numRepeats)
            {
                double x = numRepeats;
                ScaleSkew scaleSkew = msRepeatFittedScaleSkew.get(numRepeats);
                if(scaleSkew == null)
                {
                    continue;
                }
                double y = scaleSkew.scale;
                double w = Math.min(mStatsTable.getReadCount(numRepeats) / readCountSum7_15, 0.2);

                xUnweighted.add(x);
                yUnweighted.add(y);
                weights.add(w);
                skews.add(scaleSkew.skew);
            }

            if(xUnweighted.size() >= 2)
            {
                double[] y = new double[yUnweighted.size()];
                double[][] x = new double[xUnweighted.size()][2];

                for (int i = 0; i < y.length; i++) {
                    y[i] = weights.get(i) * yUnweighted.get(i);
                    x[i][0] = weights.get(i) * xUnweighted.get(i);
                    x[i][1] = weights.get(i);
                }

                OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
                regression.setNoIntercept(true);
                regression.newSampleData(y, x);

                double[] regressionParameters = regression.estimateRegressionParameters();

                // gradient must be positive
                if(Doubles.greaterThan(regressionParameters[0], 0.0))
                {
                    scaleFitGradient = regressionParameters[0];
                    scaleFitIntercept = regressionParameters[1];

                    // also calculate the microsatellite skew as the weighted average of the skews
                    microsatelliteSkew = 0.0;
                    double weightSum = 0.0;
                    for(int i = 0; i < skews.size(); ++i)
                    {
                        microsatelliteSkew += weights.get(i) * skews.get(i);
                        weightSum += weights.get(i);
                    }
                    microsatelliteSkew /= weightSum;
                }
            }
        }

        // fall back
        if(Double.isNaN(scaleFitGradient))
        {
            scaleFitGradient = relevantDefaults.slope;
            scaleFitIntercept = relevantDefaults.scale_six - scaleFitGradient * 6;
        }

        if(Double.isNaN(microsatelliteSkew))
        {
            if(Math.min(totalShortenedCount, totalLengthenedCount) < 50)
                microsatelliteSkew = 1.0;
            else
                microsatelliteSkew = Math.max(Math.log10((double)totalShortenedCount / totalLengthenedCount) + 1, 0.01);
        }

        double optimalScaleRepeat4 = relevantDefaults.scale_four;
        double optimalScaleRepeat5 = relevantDefaults.scale_five;
        double optimalScaleRepeat6 = relevantDefaults.scale_six;

        ScaleSkew scaleSkewAt4 = msRepeatFittedScaleSkew.get(4);
        if(scaleSkewAt4 != null)
        {
            optimalScaleRepeat4 = scaleSkewAt4.scale;
        }

        ScaleSkew scaleSkewAt5 = msRepeatFittedScaleSkew.get(5);
        if(scaleSkewAt5 != null)
        {
            optimalScaleRepeat5 = scaleSkewAt5.scale;
        }

        ScaleSkew scaleSkewAt6 = msRepeatFittedScaleSkew.get(6);
        if(scaleSkewAt6 != null)
        {
            optimalScaleRepeat6 = scaleSkewAt6.scale;
        }

        // assign the params
        mJitterModelParams = new JitterModelParams(
                mStatsTable.RepeatUnit, mStatsTable.ConsensusType,
                optimalScaleRepeat4, optimalScaleRepeat5, optimalScaleRepeat6,
                scaleFitGradient, scaleFitIntercept, microsatelliteSkew);
    }

    @Nullable
    public ScaleSkew gridSearch(int numRepeats, JitterTableRow row, Double fallbackScale)
    {
        double minLoss = Double.MAX_VALUE;
        double bestScale = Double.NaN;
        double bestSkew = Double.NaN;

        if(row == null)
        {
            return null;
        }

        // we need to find the fitted scale of the previous number of repeated
        double lengthMinusOneScale = 0;
        if (msRepeatFittedScaleSkew.containsKey(numRepeats - 1))
        {
            lengthMinusOneScale = msRepeatFittedScaleSkew.get(numRepeats - 1).scale;
        }

        JitterModelLoss lossCalc = new JitterModelLoss(row, numRepeats, lengthMinusOneScale);

        for(double scale = 0.05; scale <= 5.001; scale = searchIncrement(scale))
        {
            for(double skew = 0.05; skew <= 5.001; skew = searchIncrement(skew))
            {
                double loss = lossCalc.totalLoss(scale, skew);

                if(loss < minLoss)
                {
                    minLoss = loss;
                    bestScale = scale;
                    bestSkew = skew;
                }
            }
        }
        double scaleToReturn = row.totalReadCount() < 20000 && !Double.isNaN(fallbackScale) ? fallbackScale : bestScale;
        return ScaleSkew.of(scaleToReturn, bestSkew, minLoss);
    }

    private static double searchIncrement(double last)
    {
        if(Doubles.lessOrEqual(last, 0.25))
            return last + 0.01;
        if(Doubles.lessThan(last, 0.5))
            return last + 0.02;
        if(Doubles.lessOrEqual(last, 2.0))
            return last + 0.05;
        return last + 0.2;
    }
}
