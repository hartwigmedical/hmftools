package com.hartwig.hmftools.qsee.status;

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
                : QcStatus.createEmpty();
    }

    public String toString()
    {
        return String.format("%s threshold(%s%s)",
                mKey,
                mOperator != null ? mOperator.operatorString() : "",
                mThreshold
        );
    }

    public ThresholdKey key() { return mKey; }
    public ComparisonOperator operator() { return mOperator; }
    public double threshold() { return mThreshold; }
}
