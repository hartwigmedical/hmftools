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
    LOH_PERCENT             (SummaryTableGroup.GENERAL, "LOH", SourceTool.PURPLE),

    TMB_SMALL_VARIANTS      (SummaryTableGroup.TMB, "SNVs/indels per MB", SourceTool.PURPLE),
    TMB_MS_INDELS           (SummaryTableGroup.TMB, "MS indels per MB", SourceTool.PURPLE),
    TMB_STRUCTURAL_VARIANTS (SummaryTableGroup.TMB, "SVs per MB", SourceTool.PURPLE),

    TINC                    (SummaryTableGroup.CONTAMINATION, "Tumor in normal", SourceTool.PURPLE, PurpleQCStatus.TINC_FAIL_LEVEL),
    CONTAMINATION           (SummaryTableGroup.CONTAMINATION, "Other DNA", SourceTool.PURPLE, PurpleQCStatus.MAX_CONTAMINATION),

    DELETED_GENES           (SummaryTableGroup.COPY_NUMBER, "Deleted genes", SourceTool.PURPLE, (double) PurpleQCStatus.MAX_DELETED_GENES),
    UNSUPPORTED_CN_SEGMENTS (SummaryTableGroup.COPY_NUMBER, "Unsupp. segments", SourceTool.PURPLE, (double) PurpleQCStatus.MAX_UNSUPPORTED_SEGMENTS),

    MAPPED_PROPORTION       (SummaryTableGroup.MAPPING, "Mapped proportion", SourceTool.BAM_METRICS),
    MIN_COVERAGE_10         (SummaryTableGroup.MAPPING, "Coverage ≥ 10", SourceTool.BAM_METRICS),
    MIN_COVERAGE_20         (SummaryTableGroup.MAPPING, "Coverage ≥ 20", SourceTool.BAM_METRICS),
    MIN_COVERAGE_30         (SummaryTableGroup.MAPPING, "Coverage ≥ 30", SourceTool.BAM_METRICS),
    MIN_COVERAGE_60         (SummaryTableGroup.MAPPING, "Coverage ≥ 60", SourceTool.BAM_METRICS),
    MIN_COVERAGE_100        (SummaryTableGroup.MAPPING, "Coverage ≥ 100", SourceTool.BAM_METRICS),
    MIN_COVERAGE_250        (SummaryTableGroup.MAPPING, "Coverage ≥ 250", SourceTool.BAM_METRICS),
    LOW_MAP_QUAL            (SummaryTableGroup.MAPPING, "Low map quality", SourceTool.BAM_METRICS),
    LOW_BASE_QUAL           (SummaryTableGroup.MAPPING, "Low base quality", SourceTool.BAM_METRICS),
    DUPLICATE_READS         (SummaryTableGroup.MAPPING, "Duplicate reads", SourceTool.BAM_METRICS),
    DUAL_STRAND_READS       (SummaryTableGroup.MAPPING, "Dual-strand reads", SourceTool.BAM_METRICS);

    private final SummaryTableGroup mGroup;
    private final String mPlotLabel;
    private final SourceTool mSourceTool;
    private final @Nullable Double mQcThreshold;

    SummaryTableFeature(SummaryTableGroup group, String plotLabel, SourceTool sourceTool, @Nullable Double qcThreshold)
    {
        mGroup = group;
        mPlotLabel = plotLabel;
        mSourceTool = sourceTool;
        mQcThreshold = qcThreshold;
    }

    SummaryTableFeature(SummaryTableGroup group, String plotLabel, SourceTool sourceTool)
    {
        mGroup = group;
        mPlotLabel = plotLabel;
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
        String featureName = FeatureKey.formMultiFieldName("FeatureGroup", mGroup.humanReadableName(), "Key", this.toString());
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

    public String plotLabel() { return mPlotLabel; }
}
