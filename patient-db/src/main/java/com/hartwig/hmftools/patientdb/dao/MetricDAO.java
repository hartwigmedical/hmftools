package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.METRIC;

import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
import com.hartwig.hmftools.common.metrics.WGSMetrics;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class MetricDAO {

    @NotNull
    private final DSLContext context;

    MetricDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void writeMetrics(@NotNull String sample, @NotNull WGSMetricWithQC metrics) {
        deleteMetricForSample(sample);

        WGSMetrics wgsMetrics = metrics.wgsMetrics();
        Double tumorMeanCoverage = wgsMetrics.tumorMeanCoverage();
        Double tumor30xCoveragePercentage = wgsMetrics.tumor30xCoveragePercentage();
        Double tumor60xCoveragePercentage = wgsMetrics.tumor60xCoveragePercentage();

        // We only write metrics for somatic runs.
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
                METRIC.TUMORCOVERAGE60XPERCENTAGE,
                METRIC.SUFFICIENTCOVERAGE)
                .values(sample,
                        DatabaseUtil.decimal(wgsMetrics.refMeanCoverage()),
                        DatabaseUtil.decimal(wgsMetrics.ref10xCoveragePercentage()),
                        DatabaseUtil.decimal(wgsMetrics.ref20xCoveragePercentage()),
                        DatabaseUtil.decimal(tumorMeanCoverage),
                        DatabaseUtil.decimal(tumor30xCoveragePercentage),
                        DatabaseUtil.decimal(tumor60xCoveragePercentage),
                        metrics.qcMetric() ? (byte) 1 : (byte) 0)
                .execute();
    }

    void deleteMetricForSample(@NotNull String sample) {
        context.delete(METRIC).where(METRIC.SAMPLEID.eq(sample)).execute();
    }
}
