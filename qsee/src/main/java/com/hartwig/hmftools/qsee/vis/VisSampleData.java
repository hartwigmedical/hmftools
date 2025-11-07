package com.hartwig.hmftools.qsee.vis;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;

public class VisSampleData
{
    private final String mSampleId;
    private final SampleType mSampleType;
    private final Feature mFeature;
    private final double mPercentileInCohort;

    public VisSampleData(String sampleId, SampleType sampleType, Feature feature, double percentileInCohort)
    {
        mSampleId = sampleId;
        mSampleType = sampleType;
        mFeature = feature;
        mPercentileInCohort = percentileInCohort;
    }

    public String sampleId() { return mSampleId; }
    public SampleType sampleType() { return mSampleType; }
    public Feature feature() { return mFeature; }
    public double percentileInCohort() { return mPercentileInCohort; }
}
