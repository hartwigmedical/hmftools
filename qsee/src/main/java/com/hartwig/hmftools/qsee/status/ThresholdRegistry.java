package com.hartwig.hmftools.qsee.status;

import static com.hartwig.hmftools.qsee.common.SampleType.NORMAL;
import static com.hartwig.hmftools.qsee.common.SampleType.TUMOR;
import static com.hartwig.hmftools.qsee.feature.FeatureType.SUMMARY_TABLE;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.CONTAMINATION;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.DELETED_GENES;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.DUAL_STRAND_READS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.DUPLICATE_READS;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.LOW_BASE_QUAL;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.LOW_MAP_QUAL;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.MAPPED_PROPORTION;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.MEAN_COVERAGE;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_10;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_100;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_20;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_250;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_30;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.COVERAGE_ABOVE_60;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.PURITY;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.TINC;
import static com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature.UNSUPPORTED_CN_SEGMENTS;
import static com.hartwig.hmftools.qsee.status.ComparisonOperator.GREATER_THAN;
import static com.hartwig.hmftools.qsee.status.ComparisonOperator.GREATER_THAN_OR_EQUAL;
import static com.hartwig.hmftools.qsee.status.ComparisonOperator.LESS_THAN;
import static com.hartwig.hmftools.qsee.status.QcStatusType.FAIL;
import static com.hartwig.hmftools.qsee.status.QcStatusType.WARN;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.NoSuchElementException;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

public final class ThresholdRegistry
{
    private final Map<ThresholdKey, QcThreshold> mThresholds;
    private final boolean mFrozen;

    private static final double THRESHOLD_NOT_SET = Double.NaN;
    private static final double NO_THRESHOLD = Double.NaN;

    private ThresholdRegistry(Map<ThresholdKey, QcThreshold> thresholds, boolean frozen)
    {
        mThresholds = thresholds;
        mFrozen = frozen;
    }

    ThresholdRegistry()
    {
        mThresholds = new LinkedHashMap<>();
        initialise();

        mFrozen = false;
    }

    private ThresholdRegistry initialise()
    {
        setCommonThreshold(SUMMARY_TABLE, MAPPED_PROPORTION.name(), FAIL, LESS_THAN, 0.95);
        setCommonThreshold(SUMMARY_TABLE, LOW_MAP_QUAL.name(), WARN, GREATER_THAN, 0.05);
        setCommonThreshold(SUMMARY_TABLE, LOW_BASE_QUAL.name(), WARN, GREATER_THAN, 0.05);
        setCommonThreshold(SUMMARY_TABLE, DUPLICATE_READS.name(), WARN, GREATER_THAN, 0.3);
        setCommonThreshold(SUMMARY_TABLE, DUAL_STRAND_READS.name(), WARN, GREATER_THAN, 0.5);

        // PURPLE QC thresholds are not handled by Qsee
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, PURITY.name(), WARN, LESS_THAN, PurpleQCStatus.MIN_PURITY);
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, DELETED_GENES.name(), WARN, GREATER_THAN, PurpleQCStatus.MAX_DELETED_GENES);
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, UNSUPPORTED_CN_SEGMENTS.name(), WARN, GREATER_THAN, PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS);
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, CONTAMINATION.name(), FAIL, GREATER_THAN, PurpleQCStatus.MAX_CONTAMINATION);
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, TINC.name(), WARN, GREATER_THAN, PurpleQCStatus.TINC_WARN_LEVEL);
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, TINC.name(), FAIL, GREATER_THAN_OR_EQUAL, PurpleQCStatus.TINC_FAIL_LEVEL);

        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, PurpleQCStatus.FAIL_NO_TUMOR.name(), FAIL, null, NO_THRESHOLD);
        setThresholdHandledElsewhere(TUMOR, SUMMARY_TABLE, PurpleQCStatus.WARN_GENDER_MISMATCH.name(), WARN, null, NO_THRESHOLD);

        setThreshold(TUMOR, SUMMARY_TABLE, MEAN_COVERAGE.name(), WARN, LESS_THAN, 70);
        setThreshold(TUMOR, SUMMARY_TABLE, COVERAGE_ABOVE_10.name(), WARN, LESS_THAN, 0.9);
        setThreshold(TUMOR, SUMMARY_TABLE, COVERAGE_ABOVE_20.name(), WARN, LESS_THAN, 0.9);
        setThreshold(TUMOR, SUMMARY_TABLE, COVERAGE_ABOVE_30.name(), WARN, LESS_THAN, 0.9);
        setThreshold(TUMOR, SUMMARY_TABLE, COVERAGE_ABOVE_60.name(), WARN, LESS_THAN, 0.8);
        setThreshold(TUMOR, SUMMARY_TABLE, COVERAGE_ABOVE_100.name(), WARN, LESS_THAN, 0.1);
        setThreshold(TUMOR, SUMMARY_TABLE, COVERAGE_ABOVE_250.name(), WARN, LESS_THAN, THRESHOLD_NOT_SET);

        setThreshold(NORMAL, SUMMARY_TABLE, MEAN_COVERAGE.name(), WARN, LESS_THAN, 20);
        setThreshold(NORMAL, SUMMARY_TABLE, COVERAGE_ABOVE_10.name(), WARN, LESS_THAN, 0.9);
        setThreshold(NORMAL, SUMMARY_TABLE, COVERAGE_ABOVE_20.name(), WARN, LESS_THAN, 0.8);
        setThreshold(NORMAL, SUMMARY_TABLE, COVERAGE_ABOVE_30.name(), WARN, LESS_THAN, 0.4);
        setThreshold(NORMAL, SUMMARY_TABLE, COVERAGE_ABOVE_60.name(), WARN, LESS_THAN, THRESHOLD_NOT_SET);
        setThreshold(NORMAL, SUMMARY_TABLE, COVERAGE_ABOVE_100.name(), WARN, LESS_THAN, THRESHOLD_NOT_SET);
        setThreshold(NORMAL, SUMMARY_TABLE, COVERAGE_ABOVE_250.name(), WARN, LESS_THAN, THRESHOLD_NOT_SET);

        return this;
    }

    ThresholdRegistry freeze()
    {
        return new ThresholdRegistry(Collections.unmodifiableMap(mThresholds), true);
    }

    public static ThresholdRegistry createDefault()
    {
        return new ThresholdRegistry().initialise().freeze();
    }

    @VisibleForTesting
    public static ThresholdRegistry createWithoutThresholds()
    {
        ThresholdRegistry thresholdRegistry = new ThresholdRegistry().initialise();
        for(QcThreshold threshold : thresholdRegistry.mThresholds.values())
        {
            ThresholdKey key = threshold.key();
            QcThreshold noQcThreshold = QcThreshold.builder(threshold).threshold(THRESHOLD_NOT_SET).build();
            thresholdRegistry.mThresholds.put(key, noQcThreshold);
        }

        return thresholdRegistry.freeze();
    }

    private void setThreshold(
            SampleType sampleType, FeatureType featureType, String featureName, QcStatusType qcStatusType,
            ComparisonOperator operator, double thresholdValue, boolean determinedElsewhere
    ){
        QcThreshold threshold = QcThreshold.builder()
                .sampleType(sampleType)
                .featureType(featureType)
                .featureName(featureName)
                .qcStatusType(qcStatusType)
                .comparisonOperator(operator)
                .threshold(thresholdValue)
                .determinedElsewhere(determinedElsewhere)
                .build();

        mThresholds.put(threshold.key(), threshold);
    }

    private void setThreshold(
            SampleType sampleType, FeatureType featureType, String featureName, QcStatusType qcStatusType,
            ComparisonOperator operator, double thresholdValue
    ){
        setThreshold(sampleType, featureType, featureName, qcStatusType, operator, thresholdValue, false);
    }

    private void setCommonThreshold(
            FeatureType featureType, String featureName, QcStatusType qcStatusType,
            ComparisonOperator operator, double thresholdValue
    ){
        setThreshold(TUMOR, featureType, featureName, qcStatusType, operator, thresholdValue);
        setThreshold(NORMAL, featureType, featureName, qcStatusType, operator, thresholdValue);
    }

    private void setThresholdHandledElsewhere(
            SampleType sampleType, FeatureType featureType, String featureName, QcStatusType qcStatusType,
            ComparisonOperator operator, double thresholdValue)
    {
        setThreshold(sampleType, featureType, featureName, qcStatusType, operator, thresholdValue, true);
    }

    boolean containsKey(ThresholdKey key) { return mThresholds.containsKey(key); }
    void overrideThreshold(ThresholdKey key, QcThreshold threshold) { mThresholds.put(key, threshold); }

    QcThreshold getThreshold(ThresholdKey key, boolean requireFrozen)
    {
        if(requireFrozen && !mFrozen)
        {
            throw new IllegalStateException("ThresholdRegistry must be frozen before thresholds can be accessed");
        }

        if(!mThresholds.containsKey(key))
        {
            throw new NoSuchElementException(String.format("No threshold defined for %s", key));
        }

        return mThresholds.get(key);
    }

    public QcThreshold getThreshold(ThresholdKey key)
    {
        return getThreshold(key, true);
    }

    public QcThreshold getThreshold(SampleType sampleType, FeatureType featureType, String featureName, QcStatusType qcStatusType)
    {
        return getThreshold(new ThresholdKey(sampleType, featureType, featureName, qcStatusType));
    }

    public QcThreshold getThreshold(SampleType sampleType, SummaryTableFeature summaryTableFeature, QcStatusType qcStatusType)
    {
        return getThreshold(sampleType, SUMMARY_TABLE, summaryTableFeature.name(), qcStatusType);
    }
}
