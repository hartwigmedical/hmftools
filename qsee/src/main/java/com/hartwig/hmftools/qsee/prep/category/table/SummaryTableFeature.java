package com.hartwig.hmftools.qsee.prep.category.table;

import static com.hartwig.hmftools.qsee.feature.NumberFormat.LOG;
import static com.hartwig.hmftools.qsee.feature.NumberFormat.NUMBER;
import static com.hartwig.hmftools.qsee.feature.NumberFormat.PERCENT;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.qsee.feature.NumberFormat;
import com.hartwig.hmftools.qsee.feature.SourceTool;

public enum SummaryTableFeature
{
    MEAN_COVERAGE(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Mean coverage", NUMBER),
    MAPPED_PROPORTION(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Mapped proportion", PERCENT),
    COVERAGE_ABOVE_10(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 10", PERCENT),
    COVERAGE_ABOVE_20(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 20", PERCENT),
    COVERAGE_ABOVE_30(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 30", PERCENT),
    COVERAGE_ABOVE_60(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 60", PERCENT),
    COVERAGE_ABOVE_100(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 100", PERCENT),
    COVERAGE_ABOVE_250(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage > 250", PERCENT),
    LOW_BASE_QUAL(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Low qual. bases", PERCENT),
    LOW_MAP_QUAL(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Low map qual. reads", PERCENT),
    DUPLICATE_READS(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Duplicate reads", PERCENT),
    DUAL_STRAND_READS(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Dual-strand reads", PERCENT),

    PURITY(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Purity", PERCENT),
    PLOIDY(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Ploidy", NUMBER),
    LOH_PERCENT(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "LOH", PERCENT),
    DELETED_GENES(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Deleted genes", LOG),
    UNSUPPORTED_CN_SEGMENTS(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Unsupp. segments", LOG),

    TINC(SourceTool.PURPLE, SummaryTableGroup.CONTAMINATION, "Tumor in normal", PERCENT),
    CONTAMINATION(SourceTool.PURPLE, SummaryTableGroup.CONTAMINATION, "Other DNA", PERCENT),

    TMB_SMALL_VARIANTS(SourceTool.PURPLE, SummaryTableGroup.TMB, "SNVs/indels per MB", LOG),
    TMB_MS_INDELS(SourceTool.PURPLE, SummaryTableGroup.TMB, "MS indels per MB", LOG),
    TMB_STRUCTURAL_VARIANTS(SourceTool.PURPLE, SummaryTableGroup.TMB, "SVs per MB", LOG);

    private final SourceTool mSourceTool;
    private final SummaryTableGroup mGroup;
    private final String mPlotLabel;
    private final NumberFormat mNumberFormat;

    SummaryTableFeature(SourceTool sourceTool, SummaryTableGroup group, String plotLabel, NumberFormat numberFormat)
    {
        mSourceTool = sourceTool;
        mGroup = group;
        mPlotLabel = plotLabel;
        mNumberFormat = numberFormat;
    }

    public static List<SourceTool> sourceTools()
    {
        return Arrays.stream(values()).map(SummaryTableFeature::sourceTool).distinct().sorted().toList();
    }

    public SourceTool sourceTool() { return mSourceTool; }
    public SummaryTableGroup group() { return mGroup; }
    public String plotLabel() { return mPlotLabel; }
    public NumberFormat numberFormat() { return mNumberFormat; }
}
