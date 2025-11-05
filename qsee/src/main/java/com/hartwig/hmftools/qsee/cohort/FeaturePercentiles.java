package com.hartwig.hmftools.qsee.cohort;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class FeaturePercentiles
{
    private final SampleType mSampleType;
    private final FeatureKey mKey;
    private final double[] mPercentiles;
    private final double[] mRefValues;

    public FeaturePercentiles(SampleType sampleType, FeatureKey key, double[] percentiles, double[] refValues)
    {
        mSampleType = sampleType;
        mKey = key;
        mPercentiles = percentiles;
        mRefValues = refValues;
    }

    public SampleType sampleType() { return mSampleType; }
    public FeatureKey key() { return mKey; }
    public double[] percentiles() { return mPercentiles; }
    public double[] refValues() { return mRefValues; }
}
