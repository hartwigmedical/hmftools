package com.hartwig.hmftools.cup.common;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.COMBINED;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE;
import static com.hartwig.hmftools.cup.common.ClassifierType.applyMinScore;
import static com.hartwig.hmftools.cup.common.CupConstants.FEATURE_DAMPEN_FACTOR;
import static com.hartwig.hmftools.cup.common.CupConstants.MIN_CLASSIFIER_SCORE;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class CupCalcs
{
    private static final double PERCENTILE_THRESHOLD_PERC = 0.25;
    private static final double PERCENTILE_PREV_FLOOR = 0.15;
    private static final int PERCENTILE_LOW = 5;
    private static final int PERCENTILE_HIGH = 95;

    public static final double LOW_PROB_THRESHOLD = 1e-20;

    public static Map<String,Double> calcPercentilePrevalence(
            final SampleData sample, int sampleCancerCount, int cancerTypeCount,
            final Map<String,double[]> cancerPercentiles, double value,  boolean useLowThreshold)
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

        final Map<String,Double> cancerPrevalences = Maps.newHashMap();

        for(Map.Entry<String,double[]> entry : cancerPercentiles.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!checkIsValidCancerType(sample, cancerType, cancerPrevalences))
                continue;

            final double[] percentiles = entry.getValue();

            double percentile = getPercentile(percentiles, threshold);

            if(sampleCancerCount > 0 && cancerType.equals(sample.cancerType()))
            {
                // adjust percentile value to remove this sample
                if(useLowThreshold && value < threshold)
                {
                    double countBelowThreshold = percentile * sampleCancerCount;
                    percentile = max((countBelowThreshold - 1) / (sampleCancerCount - 1), 0);
                }
                else if(!useLowThreshold && value > threshold)
                {
                    double countAboveThreshold = (1 - percentile) * sampleCancerCount;
                    percentile = min(1 - (countAboveThreshold - 1) / (sampleCancerCount - 1), 1.00);
                }
            }

            percentile = useLowThreshold ? max(percentile, prevFloor) :  max(1 - percentile, prevFloor);
            cancerPrevalences.put(cancerType, percentile);
        }

        convertToPercentages(cancerPrevalences);

        return cancerPrevalences;
    }

    public static SampleResult calcCombinedFeatureResult(final SampleData sample, final List<SampleResult> allResults, boolean purgeContributors)
    {
        final Map<String,Double> cancerPrevalenceValues = Maps.newHashMap();

        final List<SampleResult> prevalenceResults = allResults.stream()
                .filter(x -> x.Result == LIKELIHOOD)
                .filter(x -> x.Category != CLASSIFIER)
                .collect(Collectors.toList());

        if(prevalenceResults.isEmpty())
            return null;

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

        dampenProbabilities(cancerPrevalenceValues, FEATURE_DAMPEN_FACTOR);
        convertToPercentages(cancerPrevalenceValues);

        // remove the contributing prevalence likelihood results
        if(purgeContributors)
            prevalenceResults.forEach(x -> allResults.remove(x));

        return new SampleResult(sample.Id, CLASSIFIER, LIKELIHOOD, FEATURE.toString(), "", cancerPrevalenceValues);
    }

    public static void convertToPercentages(final Map<String,Double> dataMap)
    {
        double valueTotal = dataMap.values().stream().mapToDouble(x -> x).sum();

        if(valueTotal == 0)
            return;

        for(Map.Entry<String,Double> entry : dataMap.entrySet())
        {
            double percentage = entry.getValue() / valueTotal;
            dataMap.put(entry.getKey(), percentage);
        }
    }

    public static void fillMissingCancerTypeValues(final Map<String,Double> dataMap, final Set<String> cancerTypes)
    {
        for(String cancerType : cancerTypes)
        {
            if(!dataMap.containsKey(cancerType))
            {
                dataMap.put(cancerType, 0.0);
            }
        }
    }

    public static void adjustLowProbabilities(final Map<String,Double> probabilityMap)
    {
        double maxProbability = probabilityMap.values().stream().mapToDouble(x -> x).max().orElse(0);

        if(maxProbability <= 0)
            return;

        // prevent values getting too small and dropping out
        if(maxProbability < LOW_PROB_THRESHOLD)
            convertToPercentages(probabilityMap);
    }

    public static void dampenProbabilities(final Map<String,Double> cancerProbabilities, double dampenFactor)
    {
        for(Map.Entry<String,Double> entry : cancerProbabilities.entrySet())
        {
            double adjustedProb = pow(entry.getValue(), dampenFactor);
            cancerProbabilities.put(entry.getKey(), adjustedProb);
        }
    }
    public static SampleResult calcCombinedClassifierScoreResult(
            final SampleData sample, final List<SampleResult> results, final String dataType, double dampenFactor)
    {
        // combined a set of classifier into a single new combined result
        final List<SampleResult> classifierResults = results.stream()
                .filter(x -> x.Category == CLASSIFIER)
                .collect(Collectors.toList());

        if(classifierResults.size() == 1)
            return null;

        final Map<String,Double> cancerTypeValues = Maps.newHashMap();

        for(final SampleResult result : classifierResults)
        {
            boolean applyMinScore = applyMinScore(ClassifierType.valueOf(result.DataType));

            for(Map.Entry<String,Double> entry : result.CancerTypeValues.entrySet())
            {
                final String cancerType = entry.getKey();

                double probability = entry.getValue();

                if(applyMinScore && sample.isCandidateCancerType(cancerType))
                    probability = max(probability, MIN_CLASSIFIER_SCORE);

                Double probTotal = cancerTypeValues.get(cancerType);

                if(probTotal == null)
                    cancerTypeValues.put(cancerType, probability);
                else
                    cancerTypeValues.put(cancerType, probTotal * probability);
            }
        }

        dampenProbabilities(cancerTypeValues, dampenFactor);
        convertToPercentages(cancerTypeValues);

        return new SampleResult(sample.Id, COMBINED, LIKELIHOOD, dataType, "", cancerTypeValues);
    }

    public static double[] adjustRefCounts(final double[] refCounts, final double[] sampleCounts, final double sampleFactor)
    {
        double[] adjustedCounts = new double[refCounts.length];

        for(int b = 0; b < refCounts.length; ++b)
        {
            adjustedCounts[b] = max(refCounts[b] - (sampleCounts[b] * sampleFactor), 0);
        }

        return adjustedCounts;
    }

}
