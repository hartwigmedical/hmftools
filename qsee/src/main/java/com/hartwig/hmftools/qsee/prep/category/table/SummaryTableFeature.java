package com.hartwig.hmftools.qsee.prep.category.table;

import static com.hartwig.hmftools.qsee.status.ComparisonOperator.GREATER_THAN;
import static com.hartwig.hmftools.qsee.status.ComparisonOperator.LESS_THAN;
import static com.hartwig.hmftools.qsee.status.QcStatusType.FAIL;
import static com.hartwig.hmftools.qsee.status.QcStatusType.WARN;

import java.util.Arrays;
import java.util.List;

import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.status.QcThreshold;

public enum SummaryTableFeature
{
    MEAN_COVERAGE           (SummaryTableGroup.GENERAL, "Mean coverage", SourceTool.BAM_METRICS),
    PURITY                  (SummaryTableGroup.GENERAL, "Purity", SourceTool.PURPLE, QcThreshold.determinedElsewhere()),
    PLOIDY                  (SummaryTableGroup.GENERAL, "Ploidy", SourceTool.PURPLE),
    LOH_PERCENT             (SummaryTableGroup.GENERAL, "LOH", SourceTool.PURPLE),

    TMB_SMALL_VARIANTS      (SummaryTableGroup.TMB, "SNVs/indels per MB", SourceTool.PURPLE),
    TMB_MS_INDELS           (SummaryTableGroup.TMB, "MS indels per MB", SourceTool.PURPLE),
    TMB_STRUCTURAL_VARIANTS (SummaryTableGroup.TMB, "SVs per MB", SourceTool.PURPLE),

    TINC                    (SummaryTableGroup.CONTAMINATION, "Tumor in normal", SourceTool.PURPLE, QcThreshold.determinedElsewhere()),
    CONTAMINATION           (SummaryTableGroup.CONTAMINATION, "Other DNA", SourceTool.PURPLE, QcThreshold.determinedElsewhere()),

    DELETED_GENES           (SummaryTableGroup.COPY_NUMBER, "Deleted genes", SourceTool.PURPLE, QcThreshold.determinedElsewhere()),
    UNSUPPORTED_CN_SEGMENTS (SummaryTableGroup.COPY_NUMBER, "Unsupp. segments", SourceTool.PURPLE, QcThreshold.determinedElsewhere()),

    MAPPED_PROPORTION       (SummaryTableGroup.MAPPING, "Mapped proportion", SourceTool.BAM_METRICS, new QcThreshold(FAIL, LESS_THAN, 0.95)),
    MIN_COVERAGE_10         (SummaryTableGroup.MAPPING, "Coverage ≥ 10", SourceTool.BAM_METRICS, new QcThreshold(WARN, LESS_THAN, 0.9)),
    MIN_COVERAGE_20         (SummaryTableGroup.MAPPING, "Coverage ≥ 20", SourceTool.BAM_METRICS, new QcThreshold(WARN, LESS_THAN, 0.9)),
    MIN_COVERAGE_30         (SummaryTableGroup.MAPPING, "Coverage ≥ 30", SourceTool.BAM_METRICS, new QcThreshold(WARN, LESS_THAN, 0.9)),
    MIN_COVERAGE_60         (SummaryTableGroup.MAPPING, "Coverage ≥ 60", SourceTool.BAM_METRICS, new QcThreshold(WARN, LESS_THAN, 0.8)),
    MIN_COVERAGE_100        (SummaryTableGroup.MAPPING, "Coverage ≥ 100", SourceTool.BAM_METRICS, new QcThreshold(WARN, LESS_THAN, 0.1)),
    MIN_COVERAGE_250        (SummaryTableGroup.MAPPING, "Coverage ≥ 250", SourceTool.BAM_METRICS),
    LOW_MAP_QUAL            (SummaryTableGroup.MAPPING, "Low map quality", SourceTool.BAM_METRICS, new QcThreshold(WARN, GREATER_THAN, 0.05)),
    LOW_BASE_QUAL           (SummaryTableGroup.MAPPING, "Low base quality", SourceTool.BAM_METRICS, new QcThreshold(WARN, GREATER_THAN, 0.05)),
    DUPLICATE_READS         (SummaryTableGroup.MAPPING, "Duplicate reads", SourceTool.BAM_METRICS, new QcThreshold(WARN, GREATER_THAN, 0.3)),
    DUAL_STRAND_READS       (SummaryTableGroup.MAPPING, "Dual-strand reads", SourceTool.BAM_METRICS, new QcThreshold(WARN, GREATER_THAN, 0.5));

    private final SummaryTableGroup mGroup;
    private final String mPlotLabel;
    private final SourceTool mSourceTool;
    private final QcThreshold mQcThreshold;

    SummaryTableFeature(SummaryTableGroup group, String plotLabel, SourceTool sourceTool, QcThreshold qcThreshold)
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

    public String plotLabel() { return mPlotLabel; }

    public boolean hasQcThreshold() { return mQcThreshold != null; }

    public QcThreshold qcThreshold()
    {
        if(mQcThreshold == null)
        {
            throw new IllegalStateException("QC threshold not set for " + this);
        }

        return mQcThreshold;
    }
}
