package com.hartwig.hmftools.patientdb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.metrics.ImmutableBamMetricSummary;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;
import com.hartwig.hmftools.patientdb.database.hmfpatients.tables.records.MetricRecord;

import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class LoadMetricsDataTest extends DatabaseTestBase
{
    @Test
    public void canWriteMetricsData()
    {
        BamMetricSummary metricSummaryTemplate = ImmutableBamMetricSummary.builder()
                .totalRegionBases(0)
                .totalReads(0)
                .duplicateReads(0)
                .dualStrandReads(0)
                .meanCoverage(0)
                .sdCoverage(0)
                .medianCoverage(0)
                .madCoverage(0)
                .lowMapQualPercent(0)
                .duplicatePercent(0)
                .unmappedPercent(0)
                .lowBaseQualPercent(0)
                .overlappingReadPercent(0)
                .cappedCoveragePercent(0)
                .build();

        BamMetricSummary tumorMetrics = metricSummaryTemplate;
        BamMetricSummary referenceMetrics = metricSummaryTemplate;

        databaseAccess.writeMetrics(TEST_SAMPLE_ID, tumorMetrics, referenceMetrics);

        List<MetricRecord> metricRecords = fetchTable(Tables.METRIC, MetricRecord.class);
        assertEquals(1, metricRecords.size());
    }
}
