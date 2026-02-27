package com.hartwig.hmftools.qsee.status;

import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

public class PurpleQCStatusConverter
{
    private final ThresholdRegistry mThresholds;

    public PurpleQCStatusConverter(ThresholdRegistry thresholds)
    {
        mThresholds = thresholds;
    }

    private QcStatus getQseeQcStatus(String featureName, QcStatusType qcStatusType)
    {
        return mThresholds.getThreshold(SampleType.TUMOR, FeatureType.SUMMARY_TABLE, featureName, qcStatusType).getQcStatus();
    }

    public QcStatus toQseeQcStatus(PurpleQCStatus purpleQCStatus)
    {
        return switch(purpleQCStatus)
        {
            case PASS -> QcStatus.createEmpty();
            case WARN_DELETED_GENES -> getQseeQcStatus(SummaryTableFeature.DELETED_GENES.name(), QcStatusType.WARN);
            case WARN_HIGH_COPY_NUMBER_NOISE -> getQseeQcStatus(SummaryTableFeature.UNSUPPORTED_CN_SEGMENTS.name(), QcStatusType.WARN);
            case WARN_LOW_PURITY -> getQseeQcStatus(SummaryTableFeature.PURITY.name(), QcStatusType.WARN);
            case WARN_TINC -> getQseeQcStatus(SummaryTableFeature.TINC.name(), QcStatusType.WARN);
            case FAIL_TINC -> getQseeQcStatus(SummaryTableFeature.TINC.name(), QcStatusType.FAIL);
            case FAIL_CONTAMINATION -> getQseeQcStatus(SummaryTableFeature.CONTAMINATION.name(), QcStatusType.FAIL);

            case FAIL_NO_TUMOR -> getQseeQcStatus(PurpleQCStatus.FAIL_NO_TUMOR.name(), QcStatusType.FAIL);
            case WARN_GENDER_MISMATCH -> getQseeQcStatus(PurpleQCStatus.WARN_GENDER_MISMATCH.name(), QcStatusType.WARN);

            default ->
                    throw new IllegalArgumentException("QC threshold not defined for PurpleQCStatus(" + purpleQCStatus + ")");
        };
    }
}
