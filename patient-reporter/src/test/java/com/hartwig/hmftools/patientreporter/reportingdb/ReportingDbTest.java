package com.hartwig.hmftools.patientreporter.reportingdb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.LimsCohort;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.ImmutableLimsCohortModel;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfigData;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortModel;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class ReportingDbTest {

    private static final String REPORTING_DB_TSV = Resources.getResource("reporting_db/reporting_db_example.tsv").getPath();

    private static final String TEST_DB_OUTPUT_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final boolean WRITE_TO_TSV = false;

    @Test
    public void canReadReportingDbTsv() throws IOException {
        List<ReportingEntry> reportingEntries = ReportingDb.read(REPORTING_DB_TSV);

        assertEquals(2, reportingEntries.size());

        ReportingEntry reportingEntry1 = reportingEntries.get(0);
        assertEquals("ABCD", reportingEntry1.tumorBarcode());
        assertEquals("sampleId1", reportingEntry1.sampleId());
        assertEquals("A", reportingEntry1.cohort());
        assertEquals("22-Oct-2019", reportingEntry1.reportDate());
        assertEquals("dna_analysis_report", reportingEntry1.reportType());
        assertEquals("0.70", reportingEntry1.purity());
        assertEquals("true", reportingEntry1.hasReliableQuality());
        assertEquals("false", reportingEntry1.hasReliablePurity());

        ReportingEntry reportingEntry2 = reportingEntries.get(1);
        assertEquals("EFGH", reportingEntry2.tumorBarcode());
        assertEquals("sampleId2", reportingEntry2.sampleId());
        assertEquals("B", reportingEntry2.cohort());
        assertEquals("23-Oct-2019", reportingEntry2.reportDate());
        assertEquals("insufficient_dna", reportingEntry2.reportType());
        assertEquals("N/A", reportingEntry2.purity());
        assertEquals("N/A", reportingEntry2.hasReliableQuality());
        assertEquals("N/A", reportingEntry2.hasReliablePurity());
    }

    @Test
    public void canWriteReportDatesToTsv() throws IOException {
        if (WRITE_TO_TSV) {
            File reportDatesTsv = new File(TEST_DB_OUTPUT_DIR + File.separator + "reporting_db_test.tsv");

            if (reportDatesTsv.createNewFile()) {
                BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTsv, true));
                writer.write("tumorBarcode\tsampleId\tcohort\treportDate\treportType\tpurity\thasReliableQuality\thasReliablePurity\n");
                writer.close();
            }

            ReportingDb.addAnalysedReportToReportingDb(reportDatesTsv.getPath(),
                    ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledInAndReliablePurity("CPCT01_SUCCESS",
                            null,
                            "CPCT"));

            ReportingDb.addQCFailReportToReportingDb(reportDatesTsv.getPath(),
                    ExampleAnalysisTestFactory.buildQCFailReport("CPCT01_FAIL",
                            QCFailReason.INSUFFICIENT_TCP_SHALLOW_WGS,
                            "CPCT"));
        }
    }

    @Test
    @Ignore
    public void canDetermineWhetherSummaryIsRequired() {
        LimsCohortModel cohortModel = buildTestCohortModel("CPCT");

        assertTrue(ReportingDb.requiresSummary(cohortModel.queryCohortData("WIDE")));
        assertFalse(ReportingDb.requiresSummary(cohortModel.queryCohortData("CPCT")));
        assertTrue(ReportingDb.requiresSummary(cohortModel.queryCohortData("CORE")));
        assertFalse(ReportingDb.requiresSummary(cohortModel.queryCohortData("CORELR02")));
        assertTrue(ReportingDb.requiresSummary(cohortModel.queryCohortData("CORELR11")));
    }

    @NotNull
    private static LimsCohortModel buildTestCohortModel(@NotNull String cohortString) {
        Map<String, LimsCohortConfigData> cohortData = Maps.newHashMap();
        LimsCohortConfigData config = ImmutableLimsCohortConfigData.builder()
                .cohortId(cohortString)
                .hospitalId(true)
                .reportGermline(false)
                .reportGermlineFlag(false)
                .reportConclusion(false)
                .reportViral(false)
                .requireHospitalId(false)
                .requireHospitalPAId(false)
                .hospitalPersonsStudy(true)
                .hospitalPersonsRequester(false)
                .outputFile(false)
                .submission(false)
                .sidePanelInfo(false)
                .build();
        cohortData.put(cohortString, config);
        return ImmutableLimsCohortModel.builder().limsCohortMap(cohortData).build();
    }
}