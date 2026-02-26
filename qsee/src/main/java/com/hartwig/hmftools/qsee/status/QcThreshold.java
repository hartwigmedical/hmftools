package com.hartwig.hmftools.qsee.status;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;

public class QcThreshold
{
    private final ThresholdKey mKey;
    private final ComparisonOperator mOperator;
    private final double mThreshold;

    public QcThreshold(ThresholdKey key, ComparisonOperator operator, double threshold)
    {
        mKey = key;
        mOperator = operator;
        mThreshold = threshold;
    }

    public QcThreshold(SampleType sampleType, FeatureType featureType, String featureName,
            QcStatusType qcStatusType, ComparisonOperator operator, double threshold)
    {
        mKey = new ThresholdKey(sampleType, featureType, featureName, qcStatusType);
        mOperator = operator;
        mThreshold = threshold;
    }

    public QcStatus getQcStatus(double sampleValue)
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
}
