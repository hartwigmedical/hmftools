package com.hartwig.hmftools.qsee.status;

import java.math.BigDecimal;

import com.hartwig.hmftools.qsee.feature.NumberFormat;

import org.jetbrains.annotations.Nullable;

public class QcStatus
{
    private final QcStatusType mType;
    @Nullable private final ComparisonOperator mOperator;
    private final double mThreshold;

    public QcStatus(QcStatusType type, @Nullable ComparisonOperator operator, double threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
    }

    public QcStatus(QcStatusType type, @Nullable ComparisonOperator operator, int threshold)
    {
        mType = type;
        mOperator = operator;
        mThreshold = threshold;
    }

    public static QcStatus createEmpty()
    {
        return new QcStatus(QcStatusType.NONE, null, Double.NaN);
    }

    public static QcStatus createPass(QcStatus qcStatus)
    {
        return new QcStatus(QcStatusType.PASS, qcStatus.operator(), qcStatus.threshold());
    }

    public QcStatusType type() { return mType; }
    @Nullable public ComparisonOperator operator() { return mOperator; }
    public double threshold() { return mThreshold; }

    @Override
    public String toString()
    {
        String operatorString = mOperator != null ? mOperator.operatorString() : "null";

        return (mType == QcStatusType.NONE)
                ? "qcStatusType(NONE)"
                : String.format("qcStatusType(%s) threshold(%s%s)", mType, operatorString, mThreshold);
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

        String operatorString = Double.isNaN(thresholdValue)
                ? ""
                : operator.operatorString();

        if(Double.isNaN(thresholdValue))
            return operatorString;

        BigDecimal value = BigDecimal.valueOf(thresholdValue);
        if(isPercent)
            value = value.multiply(BigDecimal.valueOf(100));

        value = value.stripTrailingZeros();
        String thresholdString = value.toPlainString();

        if(isPercent)
            thresholdString = thresholdString + "%";

        return operatorString + thresholdString;
    }
}
