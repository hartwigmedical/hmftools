package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.common.metrics.WGSMetricQC.hasSufficientCoverage;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.METRIC;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;

import org.jooq.DSLContext;

public class MetricDAO
{
    private final DSLContext context;

    MetricDAO(final DSLContext context)
    {
        this.context = context;
    }

    void writeMetrics(final String sample, final BamMetricSummary tumorMetrics, final BamMetricSummary refMetrics)
    {
        deleteMetricForSample(sample);

        boolean sufficientCoverage = hasSufficientCoverage(tumorMetrics, refMetrics);

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
                        DatabaseUtil.decimal(0),
                        DatabaseUtil.decimal(refMetrics.lowMapQualPercent()),
                        DatabaseUtil.decimal(refMetrics.duplicatePercent()),
                        DatabaseUtil.decimal(refMetrics.unmappedPercent()),
                        DatabaseUtil.decimal(refMetrics.lowBaseQualPercent()),
                        DatabaseUtil.decimal(refMetrics.overlappingReadPercent()),
                        DatabaseUtil.decimal(refMetrics.cappedCoveragePercent()),
                        DatabaseUtil.decimal(refMetrics.totalFilteredPercent()),
                        DatabaseUtil.decimal(refMetrics.coveragePercent(1)),
                        DatabaseUtil.decimal(refMetrics.coveragePercent(10)),
                        DatabaseUtil.decimal(refMetrics.coveragePercent(20)),
                        DatabaseUtil.decimal(refMetrics.coveragePercent(30)),
                        DatabaseUtil.decimal(refMetrics.coveragePercent(60)),
                        DatabaseUtil.decimal(tumorMetrics.meanCoverage()),
                        DatabaseUtil.decimal(tumorMetrics.sdCoverage()),
                        tumorMetrics.medianCoverage(),
                        tumorMetrics.madCoverage(),
                        DatabaseUtil.decimal(0),
                        DatabaseUtil.decimal(tumorMetrics.lowMapQualPercent()),
                        DatabaseUtil.decimal(tumorMetrics.duplicatePercent()),
                        DatabaseUtil.decimal(tumorMetrics.unmappedPercent()),
                        DatabaseUtil.decimal(tumorMetrics.lowBaseQualPercent()),
                        DatabaseUtil.decimal(tumorMetrics.overlappingReadPercent()),
                        DatabaseUtil.decimal(tumorMetrics.cappedCoveragePercent()),
                        DatabaseUtil.decimal(tumorMetrics.totalFilteredPercent()),
                        DatabaseUtil.decimal(tumorMetrics.coveragePercent(1)),
                        DatabaseUtil.decimal(tumorMetrics.coveragePercent(10)),
                        DatabaseUtil.decimal(tumorMetrics.coveragePercent(20)),
                        DatabaseUtil.decimal(tumorMetrics.coveragePercent(30)),
                        DatabaseUtil.decimal(tumorMetrics.coveragePercent(60)),
                        sufficientCoverage ? (byte) 1 : (byte) 0)
                .execute();
    }

    void deleteMetricForSample(final String sample)
    {
        context.delete(METRIC).where(METRIC.SAMPLEID.eq(sample)).execute();
    }
}
