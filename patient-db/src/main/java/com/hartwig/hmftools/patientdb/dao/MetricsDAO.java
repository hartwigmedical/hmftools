package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.METRICS;

import com.hartwig.hmftools.common.metrics.WGSMetrics;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class MetricsDAO {

    @NotNull
    private final DSLContext context;

    MetricsDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeMetrics(@NotNull String sample, @NotNull WGSMetrics metrics) {
        context.delete(METRICS).where(METRICS.SAMPLEID.eq(sample)).execute();

        Double tumorMeanCoverage = metrics.tumorMeanCoverage();
        Double tumor30xCoveragePercentage = metrics.tumor30xCoveragePercentage();
        Double tumor60xCoveragePercentage = metrics.tumor60xCoveragePercentage();

        // KODU: We only write metrics for somatic runs.
        assert tumorMeanCoverage != null;
        assert tumor30xCoveragePercentage != null;
        assert tumor60xCoveragePercentage != null;

        context.insertInto(METRICS,
                METRICS.SAMPLEID,
                METRICS.REFMEANCOVERAGE,
                METRICS.REFCOVERAGE10XPERCENTAGE,
                METRICS.REFCOVERAGE20XPERCENTAGE,
                METRICS.TUMORMEANCOVERAGE,
                METRICS.TUMORCOVERAGE30XPERCENTAGE,
                METRICS.TUMORCOVERAGE60XPERCENTAGE)
                .values(sample,
                        DatabaseUtil.decimal(metrics.refMeanCoverage()),
                        DatabaseUtil.decimal(metrics.ref10xCoveragePercentage()),
                        DatabaseUtil.decimal(metrics.ref20xCoveragePercentage()),
                        DatabaseUtil.decimal(tumorMeanCoverage),
                        DatabaseUtil.decimal(tumor30xCoveragePercentage),
                        DatabaseUtil.decimal(tumor60xCoveragePercentage))
                .execute();
    }
}
