package com.hartwig.hmftools.cobalt.norm;

import java.util.Collections;
import java.util.Map;

public class NormCalcData
{
    public final double SampleMeanReadDepth;
    public final double SampleMedianReadDepth;
    public final int SampleFilteredRegionCount;

    public final Map<Integer, Double> GcBucketMedians;

    public static NormCalcData INVALID = new NormCalcData(
            0, 0, 0, Collections.emptyMap());

    public NormCalcData(
            final double sampleMeanReadDepth, final double sampleMedianReadDepth, final int sampleFilteredRegionCount,
            final Map<Integer, Double> gcBucketMedians)
    {
        SampleMeanReadDepth = sampleMeanReadDepth;
        SampleMedianReadDepth = sampleMedianReadDepth;
        SampleFilteredRegionCount = sampleFilteredRegionCount;
        GcBucketMedians = gcBucketMedians;
    }

    public double sampleMedianNormalisation()
    {
        return SampleMeanReadDepth > 0 ? SampleMedianReadDepth / SampleMeanReadDepth : 0;
    }

}
