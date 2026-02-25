package com.hartwig.hmftools.qsee.prep.category.table;

import static com.hartwig.hmftools.qsee.feature.NumberFormat.LOG;
import static com.hartwig.hmftools.qsee.feature.NumberFormat.NUMBER;
import static com.hartwig.hmftools.qsee.feature.NumberFormat.PERCENT;

import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;

import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.NumberFormat;
import com.hartwig.hmftools.qsee.feature.PlotMetadata;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.status.QcStatus;
import com.hartwig.hmftools.qsee.status.QcThreshold;

public enum SummaryTableFeature
{
    MEAN_COVERAGE(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Mean coverage", NUMBER),
    MAPPED_PROPORTION(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Mapped proportion", PERCENT),
    COVERAGE_ABOVE_10(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage >10x", PERCENT),
    COVERAGE_ABOVE_20(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage >20x", PERCENT),
    COVERAGE_ABOVE_30(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage >30x", PERCENT),
    COVERAGE_ABOVE_60(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage >60x", PERCENT),
    COVERAGE_ABOVE_100(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage >100x", PERCENT),
    COVERAGE_ABOVE_250(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Coverage >250x", PERCENT),
    LOW_BASE_QUAL(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Low quality bases", PERCENT),
    LOW_MAP_QUAL(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Low map quality reads", PERCENT),
    DUPLICATE_READS(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Duplicate reads", PERCENT),
    DUAL_STRAND_READS(SourceTool.BAM_METRICS, SummaryTableGroup.MAPPING, "Dual-strand reads", PERCENT),

    PURITY(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Purity", PERCENT),
    PLOIDY(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Ploidy", NUMBER),
    LOH_PERCENT(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "LOH", PERCENT),
    DELETED_GENES(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Deleted genes", LOG),
    UNSUPPORTED_CN_SEGMENTS(SourceTool.PURPLE, SummaryTableGroup.COPY_NUMBER, "Unsupported CN segments", LOG),

    TINC(SourceTool.PURPLE, SummaryTableGroup.CONTAMINATION, "Tumor in normal contamination", PERCENT),
    CONTAMINATION(SourceTool.PURPLE, SummaryTableGroup.CONTAMINATION, "Other DNA contamination", PERCENT),

    TMB_SMALL_VARIANTS(SourceTool.PURPLE, SummaryTableGroup.TMB, "SNVs/indels per MB", LOG),
    TMB_MS_INDELS(SourceTool.PURPLE, SummaryTableGroup.TMB, "MS indels per MB", LOG),
    TMB_STRUCTURAL_VARIANTS(SourceTool.PURPLE, SummaryTableGroup.TMB, "SVs per MB", LOG);

    private final SourceTool mSourceTool;
    private final SummaryTableGroup mGroup;
    private final String mDisplayName;
    private final NumberFormat mNumberFormat;

    SummaryTableFeature(SourceTool sourceTool, SummaryTableGroup group, String displayName, NumberFormat numberFormat)
    {
        mSourceTool = sourceTool;
        mGroup = group;
        mDisplayName = displayName;
        mNumberFormat = numberFormat;
    }

    public static List<SourceTool> sourceTools()
    {
        return Arrays.stream(values()).map(SummaryTableFeature::sourceTool).distinct().sorted().toList();
    }

    public SourceTool sourceTool() { return mSourceTool; }
    public SummaryTableGroup group() { return mGroup; }
    public String displayName() { return mDisplayName; }
    public NumberFormat numberFormat() { return mNumberFormat; }

    public static void putFeature(EnumMap<SummaryTableFeature, Feature> featuresMap, SummaryTableFeature summaryTableFeature,
            double value, QcStatus qcStatus)
    {
        FeatureKey key = new FeatureKey(summaryTableFeature.toString(), FeatureType.SUMMARY_TABLE, summaryTableFeature.sourceTool());

        PlotMetadata metadata = PlotMetadata.builder()
                .featureGroup(summaryTableFeature.group().displayName())
                .displayName(summaryTableFeature.displayName())
                .numberFormat(summaryTableFeature.numberFormat())
                .qcStatus(qcStatus)
                .build();

        Feature feature = new Feature(key, value, metadata);
        featuresMap.put(summaryTableFeature, feature);
    }

    public static void putFeature(EnumMap<SummaryTableFeature, Feature> featuresMap, SummaryTableFeature summaryTableFeature,
            double value, QcThreshold qcThreshold)
    {
        putFeature(featuresMap, summaryTableFeature, value, qcThreshold.getQcStatus(value));
    }
}
