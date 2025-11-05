package com.hartwig.hmftools.qsee.table;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public enum SummaryTableFeature
{
    PURITY                  (SummaryTableGroup.GENERAL, "Purity", SourceTool.PURPLE),
    PLOIDY                  (SummaryTableGroup.GENERAL, "Ploidy", SourceTool.PURPLE),
    TINC                    (SummaryTableGroup.GENERAL, "TINC level", SourceTool.PURPLE),
    DELETED_GENES           (SummaryTableGroup.GENERAL, "Deleted genes", SourceTool.PURPLE),
    UNSUPPORTED_CN_SEGMENTS (SummaryTableGroup.GENERAL, "Unsupported CN segments", SourceTool.PURPLE),
    LOH_PERCENT             (SummaryTableGroup.GENERAL, "LOH percent", SourceTool.PURPLE),
    CHIMERISM_PERCENT       (SummaryTableGroup.GENERAL, "Chimerism percent", SourceTool.PURPLE),
    CONTAMINATION           (SummaryTableGroup.GENERAL, "Contamination", SourceTool.PURPLE),
    CONSANGUINITY           (SummaryTableGroup.GENERAL, "Consanguinity prop.", SourceTool.AMBER),

    TMB_SMALL_VARIANTS      (SummaryTableGroup.TMB, "Small variants per MB", SourceTool.PURPLE),
    TMB_MS_INDELS           (SummaryTableGroup.TMB, "MS indels per per MB", SourceTool.PURPLE),
    TMB_STRUCTURAL_VARIANTS (SummaryTableGroup.TMB, "Structural variants", SourceTool.PURPLE),

    MEAN_COVERAGE           (SummaryTableGroup.COVERAGE, "Mean coverage", SourceTool.BAM_METRICS),
    MIN_COVERAGE_10         (SummaryTableGroup.COVERAGE, "Coverage ≥ 10", SourceTool.BAM_METRICS),
    MIN_COVERAGE_30         (SummaryTableGroup.COVERAGE, "Coverage ≥ 30", SourceTool.BAM_METRICS),
    MIN_COVERAGE_100        (SummaryTableGroup.COVERAGE, "Coverage ≥ 100", SourceTool.BAM_METRICS),
    MIN_COVERAGE_250        (SummaryTableGroup.COVERAGE, "Coverage ≥ 250", SourceTool.BAM_METRICS),
    LOW_MAP_QUAL            (SummaryTableGroup.COVERAGE, "Low map qual percent", SourceTool.BAM_METRICS),
    LOW_BASE_QUAL           (SummaryTableGroup.COVERAGE, "Low base qual percent", SourceTool.BAM_METRICS),

    DUPLICATE_READS         (SummaryTableGroup.READ, "Duplicate reads rate", SourceTool.BAM_METRICS),
    DUAL_STRAND_READS       (SummaryTableGroup.READ, "Dual-strand reads rate", SourceTool.BAM_METRICS);

    private final SummaryTableGroup mGroup;
    private final String mMetric;
    private final SourceTool mSourceTool;

    SummaryTableFeature(SummaryTableGroup group, String metric, SourceTool sourceTool)
    {
        mGroup = group;
        mMetric = metric;
        mSourceTool = sourceTool;
    }

    public static List<SourceTool> sourceTools()
    {
        return Arrays.stream(values()).map(SummaryTableFeature::sourceTool).distinct().sorted().toList();
    }

    public SourceTool sourceTool() { return mSourceTool; }

    public FeatureKey key()
    {
        String featureName = FeatureKey.formMultiFieldName("group", mGroup.humanReadableName(), "metric", mMetric);
        return new FeatureKey(FeatureType.SUMMARY_TABLE, featureName);
    }
}
