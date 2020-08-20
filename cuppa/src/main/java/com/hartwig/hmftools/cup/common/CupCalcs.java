package com.hartwig.hmftools.cup.common;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;

import java.util.Map;

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
