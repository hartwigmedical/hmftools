package com.hartwig.hmftools.qsee.prep.category.table;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

import org.jetbrains.annotations.Nullable;

public enum SummaryTableFeature
{
    MEAN_COVERAGE           (SummaryTableGroup.GENERAL, "Mean coverage", SourceTool.BAM_METRICS),
    PURITY                  (SummaryTableGroup.GENERAL, "Purity", SourceTool.PURPLE, PurpleQCStatus.MIN_PURITY),
    PLOIDY                  (SummaryTableGroup.GENERAL, "Ploidy", SourceTool.PURPLE),
    LOH_PERCENT             (SummaryTableGroup.GENERAL, "LOH percent", SourceTool.PURPLE),

    TMB_SMALL_VARIANTS      (SummaryTableGroup.TMB, "Small variants per MB", SourceTool.PURPLE),
    TMB_MS_INDELS           (SummaryTableGroup.TMB, "Microsatellite indels per MB", SourceTool.PURPLE),
    TMB_STRUCTURAL_VARIANTS (SummaryTableGroup.TMB, "Structural variants per MB", SourceTool.PURPLE),

    TINC                    (SummaryTableGroup.CONTAMINATION, "Tumor in normal contamination", SourceTool.PURPLE, PurpleQCStatus.TINC_FAIL_LEVEL),
    CONTAMINATION           (SummaryTableGroup.CONTAMINATION, "Other DNA contamination", SourceTool.PURPLE, PurpleQCStatus.MAX_CONTAMINATION),

    DELETED_GENES           (SummaryTableGroup.COPY_NUMBER, "Deleted genes", SourceTool.PURPLE, (double) PurpleQCStatus.MAX_DELETED_GENES),
    UNSUPPORTED_CN_SEGMENTS (SummaryTableGroup.COPY_NUMBER, "Unsupported CN segments", SourceTool.PURPLE, (double) PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS),

    MAPPED_PROPORTION       (SummaryTableGroup.MAPPING, "Mapped proportion", SourceTool.BAM_METRICS),
    MIN_COVERAGE_10         (SummaryTableGroup.MAPPING, "Coverage ≥ 10", SourceTool.BAM_METRICS),
    MIN_COVERAGE_20         (SummaryTableGroup.MAPPING, "Coverage ≥ 20", SourceTool.BAM_METRICS),
    MIN_COVERAGE_30         (SummaryTableGroup.MAPPING, "Coverage ≥ 30", SourceTool.BAM_METRICS),
    MIN_COVERAGE_60         (SummaryTableGroup.MAPPING, "Coverage ≥ 60", SourceTool.BAM_METRICS),
    MIN_COVERAGE_100        (SummaryTableGroup.MAPPING, "Coverage ≥ 100", SourceTool.BAM_METRICS),
    MIN_COVERAGE_250        (SummaryTableGroup.MAPPING, "Coverage ≥ 250", SourceTool.BAM_METRICS),
    LOW_MAP_QUAL            (SummaryTableGroup.MAPPING, "Low map qual percent", SourceTool.BAM_METRICS),
    LOW_BASE_QUAL           (SummaryTableGroup.MAPPING, "Low base qual percent", SourceTool.BAM_METRICS),
    DUPLICATE_READS         (SummaryTableGroup.MAPPING, "Duplicate reads rate", SourceTool.BAM_METRICS),
    DUAL_STRAND_READS       (SummaryTableGroup.MAPPING, "Dual-strand reads rate", SourceTool.BAM_METRICS);

    private final SummaryTableGroup mGroup;
    private final String mMetric;
    private final SourceTool mSourceTool;
    private final @Nullable Double mQcThreshold;

    SummaryTableFeature(SummaryTableGroup group, String metric, SourceTool sourceTool, @Nullable Double qcThreshold)
    {
        mGroup = group;
        mMetric = metric;
        mSourceTool = sourceTool;
        mQcThreshold = qcThreshold;
    }

    SummaryTableFeature(SummaryTableGroup group, String metric, SourceTool sourceTool)
    {
        mGroup = group;
        mMetric = metric;
        mSourceTool = sourceTool;
        mQcThreshold = null;
    }

    public static List<SourceTool> sourceTools()
    {
        return Arrays.stream(values()).map(SummaryTableFeature::sourceTool).distinct().sorted().toList();
    }

    public SourceTool sourceTool() { return mSourceTool; }

    public FeatureKey key()
    {
        String featureName = FeatureKey.formMultiFieldName("FeatureGroup", mGroup.humanReadableName(), "Metric", mMetric);
        return new FeatureKey(featureName, FeatureType.SUMMARY_TABLE, mSourceTool);
    }

    public Double qcThreshold()
    {
        if(mQcThreshold == null)
        {
            throw new IllegalStateException("QC threshold not set for feature: " + this);
        }

        return mQcThreshold;
    }
}
