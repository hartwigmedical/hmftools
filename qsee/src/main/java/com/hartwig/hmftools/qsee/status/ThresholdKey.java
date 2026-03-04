package com.hartwig.hmftools.qsee.status;

import java.util.Objects;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;

import org.jetbrains.annotations.Nullable;

public class ThresholdKey
{
    private final SampleType mSampleType;
    private final FeatureType mFeatureType;
    @Nullable private final String mFeatureName;
    private final QcStatusType mQcStatusType;

    public ThresholdKey(SampleType sampleType, FeatureType featureType, @Nullable String featureName, QcStatusType qcStatusType)
    {
        mSampleType = sampleType;
        mFeatureType = featureType;
        mFeatureName = featureName;
        mQcStatusType = qcStatusType;
    }

    public SampleType sampleType(){ return mSampleType;}
    public FeatureType featureType() { return mFeatureType; }
    @Nullable public String featureName() { return mFeatureName; }
    public QcStatusType qcStatusType() { return mQcStatusType; }

    @Override
    public String toString()
    {
        return String.format("sampleType(%s) featureType(%s) featureName(%s) qcStatusType(%s)",
                mSampleType, mFeatureType, mFeatureName, mQcStatusType);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        final ThresholdKey other = (ThresholdKey) o;

        boolean featureNameEqual = Objects.equals(mFeatureName, other.mFeatureName);

        return mSampleType == other.mSampleType &&
                mFeatureType == other.mFeatureType &&
                mQcStatusType == other.mQcStatusType &&
                featureNameEqual;
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mSampleType, mFeatureType, mFeatureName, mQcStatusType);
    }
}
