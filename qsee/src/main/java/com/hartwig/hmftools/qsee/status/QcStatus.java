package com.hartwig.hmftools.qsee.status;

public class QcStatus
{
    private final QcStatusType mType;
    private final ComparisonOperator mOperator;
    private final double mThreshold;

    public QcStatus(QcStatusType type, ComparisonOperator operator, double threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
    }

    public QcStatus(QcStatusType type, ComparisonOperator operator, int threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
    }

    public static QcStatus createEmpty()
    {
        return new QcStatus(QcStatusType.NONE, null, Double.NaN);
    }

    public QcStatusType type() { return mType; }
    public ComparisonOperator operator() { return mOperator; }
    public double threshold() { return mThreshold; }

    public String toString()
    {
        return (mType == QcStatusType.NONE)
                ? "qcStatusType(NONE)"
                : String.format("qcStatusType(%s) threshold(%s%s)", mType, mOperator.operatorString(), mThreshold);
    }
}
