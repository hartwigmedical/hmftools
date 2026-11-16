package com.hartwig.hmftools.qsee.feature;

import java.util.List;

import com.hartwig.hmftools.qsee.common.SampleType;

public class SampleFeatures
{
    private final String mSampleId;
    private final SampleType mSampleType;
    private final List<Feature> mFeatures;

    public SampleFeatures(String sampleId, SampleType sampleType, List<Feature> features)
    {
        mSampleId = sampleId;
        mSampleType = sampleType;
        mFeatures = features;
    }

    public String sampleId() { return mSampleId; }
    public SampleType sampleType() { return mSampleType; }
    public List<Feature> features() { return mFeatures; }
}
