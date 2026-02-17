package com.hartwig.hmftools.qsee.feature;

public class QcStatus
{
    public final String mStatus;
    public final double mThreshold;

    public static final String NO_STATUS = "";

    public QcStatus(String status, double threshold)
    {
        mStatus = status;
        mThreshold = threshold;
    }
}
