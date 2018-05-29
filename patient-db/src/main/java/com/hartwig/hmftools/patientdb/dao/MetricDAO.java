package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.METRIC;

import com.hartwig.hmftools.common.metrics.WGSMetrics;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class MetricDAO {

    @NotNull
    private final DSLContext context;

    MetricDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeMetrics(@NotNull String sample, @NotNull WGSMetrics metrics) {
        context.delete(METRIC).where(METRIC.SAMPLEID.eq(sample)).execute();

        Double tumorMeanCoverage = metrics.tumorMeanCoverage();
        Double tumor30xCoveragePercentage = metrics.tumor30xCoveragePercentage();
        Double tumor60xCoveragePercentage = metrics.tumor60xCoveragePercentage();

        // KODU: We only write metrics for somatic runs.
        assert tumorMeanCoverage != null;
        assert tumor30xCoveragePercentage != null;
        assert tumor60xCoveragePercentage != null;

        context.insertInto(METRIC,
                METRIC.SAMPLEID,
                METRIC.REFMEANCOVERAGE,
                METRIC.REFCOVERAGE10XPERCENTAGE,
                METRIC.REFCOVERAGE20XPERCENTAGE,
                METRIC.TUMORMEANCOVERAGE,
                METRIC.TUMORCOVERAGE30XPERCENTAGE,
                METRIC.TUMORCOVERAGE60XPERCENTAGE)
                .values(sample,
                        DatabaseUtil.decimal(metrics.refMeanCoverage()),
                        DatabaseUtil.decimal(metrics.ref10xCoveragePercentage()),
                        DatabaseUtil.decimal(metrics.ref20xCoveragePercentage()),
                        DatabaseUtil.decimal(tumorMeanCoverage),
                        DatabaseUtil.decimal(tumor30xCoveragePercentage),
                        DatabaseUtil.decimal(tumor60xCoveragePercentage))
                .execute();
    }

    void deleteMetricSample(@NotNull String sample) {
        context.delete(METRIC).where(METRIC.SAMPLEID.eq(sample)).execute();
    }
}
