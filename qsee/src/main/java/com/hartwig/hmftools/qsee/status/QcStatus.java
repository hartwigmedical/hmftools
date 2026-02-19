package com.hartwig.hmftools.qsee.status;

public class QcStatus
{
    private final QcStatusType mType;
    private final ComparisonOperator mOperator;
    private final double mThreshold;

    private final boolean mThresholdIsInteger;

    public QcStatus(QcStatusType type, ComparisonOperator operator, double threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
        mThresholdIsInteger = false;
    }

    public QcStatus(QcStatusType type, ComparisonOperator operator, int threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
        mThresholdIsInteger = true;
    }

    public static QcStatus createEmpty()
    {
        return new QcStatus(QcStatusType.NONE, null, Double.NaN);
    }

    public String toString()
    {
        if(mType == QcStatusType.NONE)
            return "";

        String thresholdString = mThresholdIsInteger
                ? String.valueOf(Math.round(mThreshold))
                : String.valueOf(mThreshold);

        return String.format("%s (%s%s)", mType, mOperator.operatorString(), thresholdString);
    }
}
