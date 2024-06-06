package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.METRIC;

import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
import com.hartwig.hmftools.common.metrics.WGSMetrics;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class MetricDAO
{

    private final DSLContext context;

    MetricDAO(final DSLContext context)
    {
        this.context = context;
    }

    void writeMetrics(final String sample, final WGSMetricWithQC metrics)
    {
        deleteMetricForSample(sample);

        WGSMetrics refMetrics = metrics.refMetrics();
        WGSMetrics tumorMetrics = metrics.tumorMetrics();

        context.insertInto(METRIC,
                        METRIC.SAMPLEID,
                        METRIC.REFMEANCOVERAGE,
                        METRIC.REFSDCOVERAGE,
                        METRIC.REFMEDIANCOVERAGE,
                        METRIC.REFMADCOVERAGE,
                        METRIC.REFPCTEXCADAPTER,
                        METRIC.REFPCTEXCMAPQ,
                        METRIC.REFPCTEXCDUPE,
                        METRIC.REFPCTEXCUNPAIRED,
                        METRIC.REFPCTEXCBASEQ,
                        METRIC.REFPCTEXCOVERLAP,
                        METRIC.REFPCTEXCCAPPED,
                        METRIC.REFPCTEXCTOTAL,
                        METRIC.REFCOVERAGE1XPERCENTAGE,
                        METRIC.REFCOVERAGE10XPERCENTAGE,
                        METRIC.REFCOVERAGE20XPERCENTAGE,
                        METRIC.REFCOVERAGE30XPERCENTAGE,
                        METRIC.REFCOVERAGE60XPERCENTAGE,
                        METRIC.TUMORMEANCOVERAGE,
                        METRIC.TUMORSDCOVERAGE,
                        METRIC.TUMORMEDIANCOVERAGE,
                        METRIC.TUMORMADCOVERAGE,
                        METRIC.TUMORPCTEXCADAPTER,
                        METRIC.TUMORPCTEXCMAPQ,
                        METRIC.TUMORPCTEXCDUPE,
                        METRIC.TUMORPCTEXCUNPAIRED,
                        METRIC.TUMORPCTEXCBASEQ,
                        METRIC.TUMORPCTEXCOVERLAP,
                        METRIC.TUMORPCTEXCCAPPED,
                        METRIC.TUMORPCTEXCTOTAL,
                        METRIC.TUMORCOVERAGE1XPERCENTAGE,
                        METRIC.TUMORCOVERAGE10XPERCENTAGE,
                        METRIC.TUMORCOVERAGE20XPERCENTAGE,
                        METRIC.TUMORCOVERAGE30XPERCENTAGE,
                        METRIC.TUMORCOVERAGE60XPERCENTAGE,
                        METRIC.SUFFICIENTCOVERAGE)
                .values(sample,
                        DatabaseUtil.decimal(refMetrics.meanCoverage()),
                        DatabaseUtil.decimal(refMetrics.sdCoverage()),
                        refMetrics.medianCoverage(),
                        refMetrics.madCoverage(),
                        DatabaseUtil.decimal(refMetrics.pctExcAdapter()),
                        DatabaseUtil.decimal(refMetrics.pctExcMapQ()),
                        DatabaseUtil.decimal(refMetrics.pctExcDupe()),
                        DatabaseUtil.decimal(refMetrics.pctExcUnpaired()),
                        DatabaseUtil.decimal(refMetrics.pctExcBaseQ()),
                        DatabaseUtil.decimal(refMetrics.pctExcOverlap()),
                        DatabaseUtil.decimal(refMetrics.pctExcCapped()),
                        DatabaseUtil.decimal(refMetrics.pctExcTotal()),
                        DatabaseUtil.decimal(refMetrics.coverage1xPercentage()),
                        DatabaseUtil.decimal(refMetrics.coverage10xPercentage()),
                        DatabaseUtil.decimal(refMetrics.coverage20xPercentage()),
                        DatabaseUtil.decimal(refMetrics.coverage30xPercentage()),
                        DatabaseUtil.decimal(refMetrics.coverage60xPercentage()),
                        DatabaseUtil.decimal(tumorMetrics.meanCoverage()),
                        DatabaseUtil.decimal(tumorMetrics.sdCoverage()),
                        tumorMetrics.medianCoverage(),
                        tumorMetrics.madCoverage(),
                        DatabaseUtil.decimal(tumorMetrics.pctExcAdapter()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcMapQ()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcDupe()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcUnpaired()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcBaseQ()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcOverlap()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcCapped()),
                        DatabaseUtil.decimal(tumorMetrics.pctExcTotal()),
                        DatabaseUtil.decimal(tumorMetrics.coverage1xPercentage()),
                        DatabaseUtil.decimal(tumorMetrics.coverage10xPercentage()),
                        DatabaseUtil.decimal(tumorMetrics.coverage20xPercentage()),
                        DatabaseUtil.decimal(tumorMetrics.coverage30xPercentage()),
                        DatabaseUtil.decimal(tumorMetrics.coverage60xPercentage()),
                        metrics.sufficientCoverage() ? (byte) 1 : (byte) 0)
                .execute();
    }

    void deleteMetricForSample(final String sample)
    {
        context.delete(METRIC).where(METRIC.SAMPLEID.eq(sample)).execute();
    }
}
