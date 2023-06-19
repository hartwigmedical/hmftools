package com.hartwig.hmftools.cobalt.norm;

import java.util.Collections;
import java.util.Map;

public class NormCalcData
{
    public final double SampleMeanReadCount;
    public final double SampleMedianReadCount;
    public final int SampleFilteredRegionCount;

    public final Map<Integer,Double> GcBucketMedians;

    public static NormCalcData INVALID = new NormCalcData(
            0, 0, 0, Collections.emptyMap());

    public NormCalcData(
            final double sampleMeanReadCount, final double sampleMedianReadCount, final int sampleFilteredRegionCount,
            final Map<Integer,Double> gcBucketMedians)
    {
        SampleMeanReadCount = sampleMeanReadCount;
        SampleMedianReadCount = sampleMedianReadCount;
        SampleFilteredRegionCount = sampleFilteredRegionCount;
        GcBucketMedians = gcBucketMedians;
    }

    public double sampleMedianNormalisation() { return SampleMeanReadCount > 0 ? SampleMedianReadCount / SampleMeanReadCount : 0; }

}
