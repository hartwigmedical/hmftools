package com.hartwig.hmftools.qsee.prep.category.table;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public enum SummaryTableFeature
{
    MEAN_COVERAGE           (SummaryTableGroup.MAPPING, "Mean coverage", SourceTool.BAM_METRICS),
    MAPPED_PROPORTION       (SummaryTableGroup.MAPPING, "Mapped proportion", SourceTool.BAM_METRICS),
    COVERAGE_ABOVE_10       (SummaryTableGroup.MAPPING, "Coverage > 10", SourceTool.BAM_METRICS),
    COVERAGE_ABOVE_20       (SummaryTableGroup.MAPPING, "Coverage > 20", SourceTool.BAM_METRICS),
    COVERAGE_ABOVE_30       (SummaryTableGroup.MAPPING, "Coverage > 30", SourceTool.BAM_METRICS),
    COVERAGE_ABOVE_60       (SummaryTableGroup.MAPPING, "Coverage > 60", SourceTool.BAM_METRICS),
    COVERAGE_ABOVE_100      (SummaryTableGroup.MAPPING, "Coverage > 100", SourceTool.BAM_METRICS),
    COVERAGE_ABOVE_250      (SummaryTableGroup.MAPPING, "Coverage > 250", SourceTool.BAM_METRICS),
    LOW_BASE_QUAL           (SummaryTableGroup.MAPPING, "Low qual. bases", SourceTool.BAM_METRICS),
    LOW_MAP_QUAL            (SummaryTableGroup.MAPPING, "Low map qual. reads", SourceTool.BAM_METRICS),
    DUPLICATE_READS         (SummaryTableGroup.MAPPING, "Duplicate reads", SourceTool.BAM_METRICS),
    DUAL_STRAND_READS       (SummaryTableGroup.MAPPING, "Dual-strand reads", SourceTool.BAM_METRICS),

    PURITY                  (SummaryTableGroup.COPY_NUMBER, "Purity", SourceTool.PURPLE),
    PLOIDY                  (SummaryTableGroup.COPY_NUMBER, "Ploidy", SourceTool.PURPLE),
    LOH_PERCENT             (SummaryTableGroup.COPY_NUMBER, "LOH", SourceTool.PURPLE),
    DELETED_GENES           (SummaryTableGroup.COPY_NUMBER, "Deleted genes", SourceTool.PURPLE),
    UNSUPPORTED_CN_SEGMENTS (SummaryTableGroup.COPY_NUMBER, "Unsupp. segments", SourceTool.PURPLE),

    TINC                    (SummaryTableGroup.CONTAMINATION, "Tumor in normal", SourceTool.PURPLE),
    CONTAMINATION           (SummaryTableGroup.CONTAMINATION, "Other DNA", SourceTool.PURPLE),

    TMB_SMALL_VARIANTS      (SummaryTableGroup.TMB, "SNVs/indels per MB", SourceTool.PURPLE),
    TMB_MS_INDELS           (SummaryTableGroup.TMB, "MS indels per MB", SourceTool.PURPLE),
    TMB_STRUCTURAL_VARIANTS (SummaryTableGroup.TMB, "SVs per MB", SourceTool.PURPLE);

    private final SummaryTableGroup mGroup;
    private final String mPlotLabel;
    private final SourceTool mSourceTool;

    SummaryTableFeature(SummaryTableGroup group, String plotLabel, SourceTool sourceTool)
    {
        mGroup = group;
        mPlotLabel = plotLabel;
        mSourceTool = sourceTool;
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

    public String plotLabel() { return mPlotLabel; }
}
