package com.hartwig.hmftools.qsee.prep;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;

public class VisSampleData
{
    private final String mSampleId;
    private final SampleType mSampleType;
    private final Feature mFeature;

    public VisSampleData(String sampleId, SampleType sampleType, Feature feature)
    {
        mSampleId = sampleId;
        mSampleType = sampleType;
        mFeature = feature;
    }

    public String sampleId() { return mSampleId; }
    public SampleType sampleType() { return mSampleType; }
    public Feature feature() { return mFeature; }
}
