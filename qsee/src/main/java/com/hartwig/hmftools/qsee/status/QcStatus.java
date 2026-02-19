package com.hartwig.hmftools.qsee.status;

import com.hartwig.hmftools.common.purple.PurpleQCStatus;

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

    public static QcStatus fromPurpleQcStatus(PurpleQCStatus purpleQCStatus)
    {
        return switch(purpleQCStatus)
        {
            case PASS -> QcStatus.createEmpty();

            case WARN_DELETED_GENES ->
                new QcStatus(QcStatusType.WARN, ComparisonOperator.GREATER_THAN, PurpleQCStatus.MAX_DELETED_GENES);

            case WARN_HIGH_COPY_NUMBER_NOISE ->
                new QcStatus(QcStatusType.WARN, ComparisonOperator.GREATER_THAN, PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS);

            case WARN_LOW_PURITY ->
                new QcStatus(QcStatusType.WARN, ComparisonOperator.LESS_THAN, PurpleQCStatus.MIN_PURITY);

            case WARN_TINC ->
                new QcStatus(QcStatusType.WARN, ComparisonOperator.GREATER_THAN, PurpleQCStatus.TINC_WARN_LEVEL);

            case FAIL_TINC ->
                new QcStatus(QcStatusType.FAIL, ComparisonOperator.GREATER_THAN_OR_EQUAL, PurpleQCStatus.TINC_FAIL_LEVEL);

            case FAIL_CONTAMINATION ->
                new QcStatus(QcStatusType.FAIL, ComparisonOperator.GREATER_THAN, PurpleQCStatus.MAX_CONTAMINATION);

            case FAIL_NO_TUMOR -> new QcStatus(QcStatusType.FAIL, null, Double.NaN);

            case WARN_GENDER_MISMATCH -> new QcStatus(QcStatusType.WARN, null, Double.NaN);

            default ->
                throw new IllegalArgumentException("QC threshold not defined for PurpleQCStatus(" + purpleQCStatus + ")");
        };
    }
}
