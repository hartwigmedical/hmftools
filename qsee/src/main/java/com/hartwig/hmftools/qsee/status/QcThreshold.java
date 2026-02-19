package com.hartwig.hmftools.qsee.status;

public class QcThreshold
{
    private final QcStatusType mType;
    private final ComparisonOperator mOperator;
    private final double mThreshold;

    public QcThreshold(QcStatusType type, ComparisonOperator operator, double threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
    }

    public QcStatus getQcStatus(double sampleValue)
    {
        boolean sampleFailsThreshold = switch(mOperator)
        {
            case LESS_THAN -> sampleValue < mThreshold;
            case LESS_THAN_OR_EQUAL -> sampleValue <= mThreshold;
            case GREATER_THAN_OR_EQUAL -> sampleValue >= mThreshold;
            case GREATER_THAN -> sampleValue > mThreshold;
        };

        return sampleFailsThreshold
                ? new QcStatus(mType, mOperator, mThreshold)
                : QcStatus.createEmpty();
    }

    public static QcThreshold determinedElsewhere()
    {
        return null;
    }

    public String toString()
    {
        return String.format("type(%s) operator(%s) threshold(%s)", mType, mOperator.operatorString(), mThreshold);
    }
}
