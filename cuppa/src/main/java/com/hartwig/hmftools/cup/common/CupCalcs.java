package com.hartwig.hmftools.cup.common;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE_PREVALENCE;
import static com.hartwig.hmftools.cup.common.ClassifierType.PERCENTILES;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class CupCalcs
{
    private static double CONTRIB_CONSTANT = 0.1;
    private static double CONTRIB_EXPONENT = 0.3;
    private static double MIN_CONTRIB = 0.02;

    public static double calcPercentileWeight(double percentile)
    {
        // PercentileContribution=MAX(IF(AND(Percentile>=0, Percentile <=1),contributionConstant+(1-2*ABS(0.5-Percentile))^contributionExponent,
        // contributionConstant/ABS(Percentile)),minContribution)

        if(percentile < 0 || percentile > 1)
        {
            return max(CONTRIB_CONSTANT / abs(percentile), MIN_CONTRIB);
        }

        return max(CONTRIB_CONSTANT + pow(1 - 2 * abs(0.5 - percentile), CONTRIB_EXPONENT), CONTRIB_CONSTANT);

        // Prediction = PRODUCT(PercenitleContribution) / SUM[PRODUCT(PercenitleContribution)]
    }

    public static void addPercentileClassifier(final SampleData sample, final List<SampleResult> results)
    {
        if(results.isEmpty())
            return;

        final Map<String,Double> cancerPercentileCalcs = Maps.newHashMap();

        for(final SampleResult result : results)
        {
            if(result.ResultType != PERCENTILE)
                continue;

            for(Map.Entry<String,Double> entry : result.CancerTypeValues.entrySet())
            {
                final String cancerType = entry.getKey();

                if(!sample.isCandidateCancerType(cancerType))
                {
                    cancerPercentileCalcs.put(cancerType, 0.0);
                    continue;
                }

                double percentile = entry.getValue();
                double percentileWeight = calcPercentileWeight(percentile);

                Double combinedWeight = cancerPercentileCalcs.get(cancerType);

                if(combinedWeight == null)
                    cancerPercentileCalcs.put(cancerType, percentileWeight);
                else
                    cancerPercentileCalcs.put(cancerType, percentileWeight * combinedWeight);
            }
        }

        double combinedWeightTotal = cancerPercentileCalcs.values().stream().mapToDouble(x -> x).sum();
        final Map<String,Double> cancerProbabilities = Maps.newHashMap();

        for(Map.Entry<String,Double> entry : cancerPercentileCalcs.entrySet())
        {
            cancerProbabilities.put(entry.getKey(), entry.getValue()/combinedWeightTotal);
        }

        SampleResult result = new SampleResult(
                results.get(0).SampleId, CLASSIFIER, LIKELIHOOD,
                ClassifierType.displayString(PERCENTILES), String.format("%.4f", combinedWeightTotal), cancerProbabilities);

        results.add(result);
    }

    private static final double PERCENTILE_THRESHOLD_PERC = 0.25;
    private static final double PERCENTILE_PREV_FLOOR = 0.15;
    private static final int PERCENTILE_LOW = 5;
    private static final int PERCENTILE_HIGH = 95;

    public static Map<String,Double> calcPercentilePrevalence(
            final Map<String,double[]> cancerPercentiles, double value, int cancerTypeCount, boolean useLowThreshold)
    {
        /*
        Low Test = % of samples below the max(observed value + 25%,min 5th percentile value in any cancerType)
        High Test = % of samples above the min(observed value - 33%,max 95th percentile value in any cancer type)
        Floor each % at 0.15/ # tumor types as for drivers
         */

        double threshold = useLowThreshold ?
                cancerPercentiles.values().stream().mapToDouble(x -> x[PERCENTILE_LOW]).min().orElse(0) :
                cancerPercentiles.values().stream().mapToDouble(x -> x[PERCENTILE_HIGH]).max().orElse(0);

        double prevFloor = PERCENTILE_PREV_FLOOR / cancerTypeCount;

        threshold = useLowThreshold ?
                max((1 + PERCENTILE_THRESHOLD_PERC) * value, threshold) : min((1 - PERCENTILE_THRESHOLD_PERC) * value, threshold);

        double prevTotal = 0;
        final Map<String,Double> cancerPrevalences = Maps.newHashMap();

        for(Map.Entry<String,double[]> entry : cancerPercentiles.entrySet())
        {
            final String cancerType = entry.getKey();
            final double[] percentiles = entry.getValue();

            double percentile = getPercentile(percentiles, threshold);

            percentile = useLowThreshold ? max(percentile, prevFloor) :  max(1 - percentile, prevFloor);

            prevTotal += percentile;
            cancerPrevalences.put(cancerType, percentile);
        }

        for(Map.Entry<String,Double> entry : cancerPrevalences.entrySet())
        {
            cancerPrevalences.put(entry.getKey(), entry.getValue() / prevTotal);
        }

        return cancerPrevalences;
    }

}
