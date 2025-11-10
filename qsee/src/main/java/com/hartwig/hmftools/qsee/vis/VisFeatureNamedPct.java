package com.hartwig.hmftools.qsee.vis;

import java.util.EnumMap;
import java.util.Map;

import com.hartwig.hmftools.qsee.cohort.FeaturePercentiles;
import com.hartwig.hmftools.qsee.cohort.NamedPercentile;
import com.hartwig.hmftools.qsee.cohort.PercentileTransformer;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureKey;

public class VisFeatureNamedPct
{
    private final SampleType mSampleType;
    private final FeatureKey mFeatureKey;
    private final Map<NamedPercentile, Double> mPctRefValueMap = new EnumMap<>(NamedPercentile.class);

    public VisFeatureNamedPct(SampleType sampleType, FeatureKey featureKey, FeaturePercentiles featurePercentiles)
    {
        mSampleType = sampleType;
        mFeatureKey = featureKey;

        for(NamedPercentile namedPercentile : NamedPercentile.values())
        {
            PercentileTransformer transformer = featurePercentiles.transformer();
            double featureValue = transformer.percentileToFeatureValue(namedPercentile.percentile());
            mPctRefValueMap.put(namedPercentile, featureValue);
        }
    }

    public SampleType sampleType() { return mSampleType; }
    public FeatureKey featureKey() { return mFeatureKey; }
    public Map<NamedPercentile, Double> pctRefValues() { return mPctRefValueMap; }
}
