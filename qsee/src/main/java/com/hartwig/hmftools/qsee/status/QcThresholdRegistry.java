package com.hartwig.hmftools.qsee.status;

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
import static com.hartwig.hmftools.qsee.status.ComparisonOperator.LESS_THAN;
import static com.hartwig.hmftools.qsee.status.QcStatusType.FAIL;
import static com.hartwig.hmftools.qsee.status.QcStatusType.WARN;

import java.util.Collections;
import java.util.EnumMap;
import java.util.Map;
import java.util.NoSuchElementException;

import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.prep.category.table.SummaryTableFeature;

public final class QcThresholdRegistry
{
    // TODO: Add thresholds for targeted mode

    private static final Map<SummaryTableFeature, QcThreshold> THRESHOLDS_COMMON;
    static
    {
        Map<SummaryTableFeature, QcThreshold> thresholds = new EnumMap<>(SummaryTableFeature.class);

        thresholds.put(MAPPED_PROPORTION, new QcThreshold(FAIL, LESS_THAN, 0.95));
        thresholds.put(LOW_MAP_QUAL, new QcThreshold(WARN, GREATER_THAN, 0.05));
        thresholds.put(LOW_BASE_QUAL, new QcThreshold(WARN, GREATER_THAN, 0.05));
        thresholds.put(DUPLICATE_READS, new QcThreshold(WARN, GREATER_THAN, 0.3));
        thresholds.put(DUAL_STRAND_READS, new QcThreshold(WARN, GREATER_THAN, 0.5));

        // Determined by PURPLE
        thresholds.put(PURITY, QcThreshold.determinedElsewhere());
        thresholds.put(TINC, QcThreshold.determinedElsewhere());
        thresholds.put(DELETED_GENES, QcThreshold.determinedElsewhere());
        thresholds.put(UNSUPPORTED_CN_SEGMENTS, QcThreshold.determinedElsewhere());
        thresholds.put(CONTAMINATION, QcThreshold.determinedElsewhere());

        THRESHOLDS_COMMON = Collections.unmodifiableMap(thresholds);
    }

    private static final Map<SummaryTableFeature, QcThreshold> THRESHOLDS_TUMOR;
    static
    {
        Map<SummaryTableFeature, QcThreshold> thresholds = new EnumMap<>(SummaryTableFeature.class);

        thresholds.put(MEAN_COVERAGE, new QcThreshold(WARN, LESS_THAN, 70));
        thresholds.put(COVERAGE_ABOVE_10, new QcThreshold(WARN, LESS_THAN, 0.9));
        thresholds.put(COVERAGE_ABOVE_20, new QcThreshold(WARN, LESS_THAN, 0.9));
        thresholds.put(COVERAGE_ABOVE_30, new QcThreshold(WARN, LESS_THAN, 0.9));
        thresholds.put(COVERAGE_ABOVE_60, new QcThreshold(WARN, LESS_THAN, 0.8));
        thresholds.put(COVERAGE_ABOVE_100, new QcThreshold(WARN, LESS_THAN, 0.1));
        thresholds.put(COVERAGE_ABOVE_250, QcThreshold.notSet());

        THRESHOLDS_TUMOR = Collections.unmodifiableMap(thresholds);
    }

    private static final Map<SummaryTableFeature, QcThreshold> THRESHOLDS_NORMAL;
    static
    {
        Map<SummaryTableFeature, QcThreshold> thresholds = new EnumMap<>(SummaryTableFeature.class);

        thresholds.put(MEAN_COVERAGE, new QcThreshold(WARN, LESS_THAN, 20));
        thresholds.put(COVERAGE_ABOVE_10, new QcThreshold(WARN, LESS_THAN, 0.9));
        thresholds.put(COVERAGE_ABOVE_20, new QcThreshold(WARN, LESS_THAN, 0.8));
        thresholds.put(COVERAGE_ABOVE_30, new QcThreshold(WARN, LESS_THAN, 0.4));
        thresholds.put(COVERAGE_ABOVE_60, QcThreshold.notSet());
        thresholds.put(COVERAGE_ABOVE_100, QcThreshold.notSet());
        thresholds.put(COVERAGE_ABOVE_250, QcThreshold.notSet());

        THRESHOLDS_NORMAL = Collections.unmodifiableMap(thresholds);
    }

    public static QcThreshold getThreshold(ThresholdGroup thresholdGroup, SummaryTableFeature summaryTableFeature)
    {
        Map<SummaryTableFeature, QcThreshold> thresholds = switch(thresholdGroup)
        {
            case COMMON -> THRESHOLDS_COMMON;
            case TUMOR -> THRESHOLDS_TUMOR;
            case NORMAL -> THRESHOLDS_NORMAL;
        };

        if(!thresholds.containsKey(summaryTableFeature))
        {
            throw new NoSuchElementException(String.format("Threshold not defined for thresholdGroup(%s) feature(%s)",
                    thresholdGroup, summaryTableFeature));
        }

        return thresholds.get(summaryTableFeature);
    }

    public static QcThreshold getCommonThreshold(SummaryTableFeature summaryTableFeature)
    {
        return getThreshold(ThresholdGroup.COMMON, summaryTableFeature);
    }

    public static QcThreshold getThreshold(SampleType sampleType, SummaryTableFeature summaryTableFeature)
    {
        ThresholdGroup thresholdGroup = switch(sampleType)
        {
            case TUMOR -> ThresholdGroup.TUMOR;
            case NORMAL -> ThresholdGroup.NORMAL;
        };

        return getThreshold(thresholdGroup, summaryTableFeature);
    }

}
