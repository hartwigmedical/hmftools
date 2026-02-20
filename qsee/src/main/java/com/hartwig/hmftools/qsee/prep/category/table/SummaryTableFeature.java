package com.hartwig.hmftools.qsee.prep.category.table;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public enum SummaryTableFeature
{
    MEAN_COVERAGE           (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Mean coverage"),
    MAPPED_PROPORTION       (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Mapped proportion"),
    COVERAGE_ABOVE_10       (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 10"),
    COVERAGE_ABOVE_20       (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 20"),
    COVERAGE_ABOVE_30       (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 30"),
    COVERAGE_ABOVE_60       (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 60"),
    COVERAGE_ABOVE_100      (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 100"),
    COVERAGE_ABOVE_250      (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 250"),
    LOW_BASE_QUAL           (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Low qual. bases"),
    LOW_MAP_QUAL            (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Low map qual. reads"),
    DUPLICATE_READS         (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Duplicate reads"),
    DUAL_STRAND_READS       (SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Dual-strand reads"),

    PURITY                  (SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Purity"),
    PLOIDY                  (SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Ploidy"),
    LOH_PERCENT             (SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "LOH"),
    DELETED_GENES           (SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Deleted genes"),
    UNSUPPORTED_CN_SEGMENTS (SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Unsupp. segments"),

    TINC                    (SourceTool.PURPLE, SummaryTableGroup.CONTAMINATION, "Tumor in normal"),
    CONTAMINATION           (SourceTool.PURPLE, SummaryTableGroup.CONTAMINATION, "Other DNA"),

    TMB_SMALL_VARIANTS      (SourceTool.PURPLE, SummaryTableGroup.TMB, "SNVs/indels per MB"),
    TMB_MS_INDELS           (SourceTool.PURPLE, SummaryTableGroup.TMB, "MS indels per MB"),
    TMB_STRUCTURAL_VARIANTS (SourceTool.PURPLE, SummaryTableGroup.TMB, "SVs per MB");

    private final SourceTool mSourceTool;
    private final SummaryTableGroup mGroup;
    private final String mPlotLabel;

    SummaryTableFeature(SourceTool sourceTool, SummaryTableGroup group, String plotLabel)
    {
        mSourceTool = sourceTool;
        mGroup = group;
        mPlotLabel = plotLabel;
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
