package com.hartwig.hmftools.qsee.cohort;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class FeaturePercentiles
{
    private final SampleType mSampleType;
    private final FeatureKey mFeatureKey;
    private final double[] mPercentiles;
    private final double[] mRefValues;

    private PercentileTransformer mTransformer = null;

    public FeaturePercentiles(SampleType sampleType, FeatureKey featureKey, double[] percentiles, double[] refValues)
    {
        mSampleType = sampleType;
        mFeatureKey = featureKey;
        mPercentiles = percentiles;
        mRefValues = refValues;
    }

    public SampleType sampleType() { return mSampleType; }
    public FeatureKey featureKey() { return mFeatureKey; }
    public double[] percentiles() { return mPercentiles; }
    public double[] refValues() { return mRefValues; }

    public PercentileTransformer transformer()
    {
        if(mTransformer == null)
            mTransformer = PercentileTransformer.fromPrefitData(mPercentiles, mRefValues);

        return mTransformer;
    }
}
