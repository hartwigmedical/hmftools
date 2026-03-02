package com.hartwig.hmftools.qsee.status;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QcThreshold
{
    private final ThresholdKey mKey;
    @Nullable private final ComparisonOperator mOperator;
    private final double mThreshold;

    private final boolean mDeterminedElsewhere;

    private QcThreshold(ThresholdKey key, @Nullable ComparisonOperator operator, double threshold, boolean determinedElsewhere)
    {
        mKey = key;
        mOperator = operator;
        mThreshold = threshold;
        mDeterminedElsewhere = determinedElsewhere;
    }

    public static Builder builder() { return new Builder(); }

    public static Builder builder(ThresholdKey key)
    {
        return builder()
                .sampleType(key.sampleType())
                .featureType(key.featureType())
                .featureName(key.featureName())
                .qcStatusType(key.qcStatusType());
    }

    public static Builder builder(QcThreshold threshold)
    {
        return builder()
                .sampleType(threshold.key().sampleType())
                .featureType(threshold.key().featureType())
                .featureName(threshold.key().featureName())
                .qcStatusType(threshold.key().qcStatusType())
                .comparisonOperator(threshold.operator())
                .threshold(threshold.threshold())
                .determinedElsewhere(threshold.determinedElsewhere());
    }

    public QcStatus getQcStatus()
    {
        return new QcStatus(mKey.qcStatusType(), mOperator, mThreshold);
    }

    public QcStatus getSampleQcStatus(double sampleValue)
    {
        if(mKey.qcStatusType() == QcStatusType.NONE)
        {
            return QcStatus.createEmpty();
        }

        boolean sampleFailsThreshold = switch(mOperator)
        {
            case LESS_THAN -> sampleValue < mThreshold;
            case LESS_THAN_OR_EQUAL -> sampleValue <= mThreshold;
            case GREATER_THAN_OR_EQUAL -> sampleValue >= mThreshold;
            case GREATER_THAN -> sampleValue > mThreshold;
        };

        return sampleFailsThreshold
                ? new QcStatus(mKey.qcStatusType(), mOperator, mThreshold)
                : new QcStatus(QcStatusType.PASS, mOperator, mThreshold);
    }

    @Override
    public String toString()
    {
        return String.format("%s threshold(%s %s)",
                mKey,
                mOperator != null ? mOperator.operatorString() : "",
                mThreshold
        );
    }

    public ThresholdKey key() { return mKey; }
    public ComparisonOperator operator() { return mOperator; }
    public double threshold() { return mThreshold; }
    public boolean determinedElsewhere() { return mDeterminedElsewhere; }

    @Override
    public boolean equals(Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        QcThreshold other = (QcThreshold) o;

        boolean thresholdValuesEqual = mThreshold == other.mThreshold || Double.isNaN(mThreshold) && Double.isNaN(other.mThreshold);

        return mKey.equals(other.mKey) &&
                mOperator == other.mOperator &&
                thresholdValuesEqual;
    }

    @Override
    public int hashCode()
    {
        return mKey.hashCode();
    }

    public static class Builder
    {
        @NotNull private SampleType mSampleType;

        @NotNull private FeatureType mFeatureType;
        @Nullable private String mFeatureName;
        @NotNull private QcStatusType mQcStatusType;

        @Nullable private ComparisonOperator mOperator;
        @NotNull private double mThreshold;

        @NotNull private boolean mDeterminedElsewhere;

        public Builder sampleType(SampleType sampleType)
        {
            mSampleType = sampleType;
            return this;
        }

        public Builder featureType(FeatureType featureType)
        {
            mFeatureType = featureType;
            return this;
        }

        public Builder featureName(@Nullable String featureName)
        {
            mFeatureName = featureName;
            return this;
        }

        public Builder qcStatusType(QcStatusType qcStatusType)
        {
            mQcStatusType = qcStatusType;
            return this;
        }

        public Builder comparisonOperator(@Nullable ComparisonOperator operator)
        {
            mOperator = operator;
            return this;
        }

        public Builder threshold(double threshold)
        {
            mThreshold = threshold;
            return this;
        }

        public Builder determinedElsewhere(boolean determinedElsewhere)
        {
            mDeterminedElsewhere = determinedElsewhere;
            return this;
        }

        public QcThreshold build()
        {
            ThresholdKey key = new ThresholdKey(mSampleType, mFeatureType, mFeatureName, mQcStatusType);
            return new QcThreshold(key, mOperator, mThreshold, mDeterminedElsewhere);
        }
    }
}
