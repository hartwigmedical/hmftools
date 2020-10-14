package com.hartwig.hmftools.cup.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ClassifierType.COMBINED;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE_PREVALENCE;
import static com.hartwig.hmftools.cup.common.CupConstants.MIN_CLASSIFIER_SCORE;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class CupCalcs
{
    private static final double PERCENTILE_THRESHOLD_PERC = 0.25;
    private static final double PERCENTILE_PREV_FLOOR = 0.15;
    private static final int PERCENTILE_LOW = 5;
    private static final int PERCENTILE_HIGH = 95;

    public static Map<String,Double> calcPercentilePrevalence(
            final Map<String,double[]> cancerPercentiles, double value, int cancerTypeCount, boolean useLowThreshold)
    {
        /*
        Low Test = % of samples below the max(observed value + 25%, min 5th percentile value in any cancerType)
        High Test = % of samples above the min(observed value - 25%, max 95th percentile value in any cancer type)
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

    public static SampleResult calcCombinedFeatureResult(final SampleData sample, final List<SampleResult> allResults)
    {
        final Map<String,Double> cancerPrevalenceValues = Maps.newHashMap();

        final List<SampleResult> prevalenceResults = allResults.stream()
                .filter(x -> x.ResultType == LIKELIHOOD)
                .filter(x -> x.Category != CLASSIFIER)
                .collect(Collectors.toList());

        for(final SampleResult result : prevalenceResults)
        {
            for(Map.Entry<String,Double> entry : result.CancerTypeValues.entrySet())
            {
                final String cancerType = entry.getKey();
                double prevalence = entry.getValue();

                Double prevTotal = cancerPrevalenceValues.get(cancerType);

                if(prevTotal == null)
                    cancerPrevalenceValues.put(cancerType, prevalence);
                else
                    cancerPrevalenceValues.put(cancerType, prevTotal * prevalence);
            }
        }

        double totalPrevalence = cancerPrevalenceValues.values().stream().mapToDouble(x -> x).sum();

        if(totalPrevalence == 0)
            return null;

        final Map<String,Double> cancerTypeValues = Maps.newHashMap();

        for(Map.Entry<String,Double> entry : cancerPrevalenceValues.entrySet())
        {
            double probability = entry.getValue() / totalPrevalence;
            cancerTypeValues.put(entry.getKey(), probability);
        }

        // remove the contributing prevalence likelihood results
        prevalenceResults.forEach(x -> allResults.remove(x));

        return new SampleResult(sample.Id, CLASSIFIER, LIKELIHOOD, FEATURE_PREVALENCE.toString(), "", cancerTypeValues);
    }

    public static SampleResult calcClassifierScoreResult(final SampleData sample, final List<SampleResult> results)
    {
        final List<SampleResult> classifierResults = results.stream()
                .filter(x -> x.Category == CLASSIFIER)
                .collect(Collectors.toList());

        final Map<String,Double> cancerTypeValues = Maps.newHashMap();

        for(final SampleResult result : classifierResults)
        {
            for(Map.Entry<String,Double> entry : result.CancerTypeValues.entrySet())
            {
                final String cancerType = entry.getKey();
                double probability = max(entry.getValue(), MIN_CLASSIFIER_SCORE);

                Double probTotal = cancerTypeValues.get(cancerType);

                if(probTotal == null)
                    cancerTypeValues.put(cancerType, probability);
                else
                    cancerTypeValues.put(cancerType, probTotal * probability);
            }
        }

        double totalProbability = cancerTypeValues.values().stream().mapToDouble(x -> x).sum();

        if(totalProbability == 0)
            return null;

        for(Map.Entry<String,Double> entry : cancerTypeValues.entrySet())
        {
            double probability = entry.getValue() / totalProbability;
            cancerTypeValues.put(entry.getKey(), probability);
        }

        return new SampleResult(sample.Id, CLASSIFIER, LIKELIHOOD, COMBINED.toString(), "", cancerTypeValues);
    }

}
