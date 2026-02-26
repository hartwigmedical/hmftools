package com.hartwig.hmftools.qsee.status;

import com.hartwig.hmftools.qsee.feature.NumberFormat;

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

    @Override
    public String toString()
    {
        return (mType == QcStatusType.NONE)
                ? "qcStatusType(NONE)"
                : String.format("qcStatusType(%s) threshold(%s%s)", mType, mOperator.operatorString(), mThreshold);
    }

    public String displayString(NumberFormat numberFormat)
    {
        return formDiplayString(mType, mOperator, mThreshold, numberFormat);
    }

    public static String formDiplayString(
            QcStatusType qcStatusType,
            ComparisonOperator operator,
            double thresholdValue,
            NumberFormat numberFormat)
    {
        if(qcStatusType == QcStatusType.NONE)
            return "";

        boolean isPercent = numberFormat == NumberFormat.PERCENT;
        if(isPercent)
            thresholdValue = thresholdValue * 100;

        boolean isInteger = thresholdValue % 1 == 0;
        String thresholdString = isInteger
                ? String.valueOf((int) thresholdValue)
                : String.valueOf(thresholdValue);

        if(isPercent && !Double.isNaN(thresholdValue))
            thresholdString = thresholdString + "%";

        String operatorString = Double.isNaN(thresholdValue)
                ? ""
                : operator.operatorString();

        return operatorString + thresholdString;
    }
}
