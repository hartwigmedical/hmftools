package com.hartwig.hmftools.patientreporter.reportingdb;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.LimsCoreCohort;
import com.hartwig.hmftools.common.lims.LimsStudy;
import com.hartwig.hmftools.patientreporter.ExampleAnalysisTestFactory;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.junit.Test;

public class ReportingDbTest {

    private static final String REPORT_DATES_TSV = Resources.getResource("reporting_db/reporting_db_example.tsv").getPath();

    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final boolean WRITE_TO_TSV = false;

    @Test
    public void canReadReportDatesTsv() throws IOException {
        List<ReportingEntry> reportingEntries = ReportingDb.read(REPORT_DATES_TSV);

        assertEquals(2, reportingEntries.size());

        ReportingEntry reportingEntry1 = reportingEntries.get(0);
        assertEquals("ABCD", reportingEntry1.tumorBarcode());
        assertEquals("sampleId1", reportingEntry1.sampleId());
        assertEquals("22-Oct-2019", reportingEntry1.reportDate());
        assertEquals("sequence_report", reportingEntry1.reportType());
        assertEquals("0.70", reportingEntry1.purity());
        assertEquals("true", reportingEntry1.hasReliableQuality());
        assertEquals("false", reportingEntry1.hasReliablePurity());
        assertEquals("A", reportingEntry1.cohort());

        ReportingEntry reportingEntry2 = reportingEntries.get(1);
        assertEquals("EFGH", reportingEntry2.tumorBarcode());
        assertEquals("sampleId2", reportingEntry2.sampleId());
        assertEquals("22-Oct-2019", reportingEntry2.reportDate());
        assertEquals("shallow_seq_low_purity", reportingEntry2.reportType());
        assertEquals("N/A", reportingEntry2.purity());
        assertEquals("N/A", reportingEntry2.hasReliableQuality());
        assertEquals("N/A", reportingEntry2.hasReliablePurity());
        assertEquals("B", reportingEntry2.cohort());

    }

    @Test
    public void canWriteReportDatesToTsv() throws IOException {
        if (WRITE_TO_TSV) {
            File reportDatesTsv = new File(REPORT_BASE_DIR + File.separator + "reporting_db_test.tsv");

            if (reportDatesTsv.createNewFile()) {
                BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTsv, true));
                writer.write("tumorBarcode\tsampleId\treportDate\treportType\tpurity\thasReliableQuality\thasReliablePurity\tcohort\n");
                writer.close();
            }

            ReportingDb.addSequenceReportToReportingDb(reportDatesTsv.getPath(), ExampleAnalysisTestFactory.buildCOLO829());

            ReportingDb.addQCFailReportToReportingDb(reportDatesTsv.getPath(),
                    ExampleAnalysisTestFactory.buildQCFailReport("LowMolecularTumorSample", QCFailReason.SHALLOW_SEQ_LOW_PURITY));
        }
    }

    @Test
    public void determineForCheckSummary() {
        String sampleId1 = "WIDE01000001T";
        LimsCoreCohort coreCohort1 = LimsCoreCohort.fromSampleId(sampleId1);
        LimsStudy study1 = LimsStudy.fromSampleId(sampleId1);
        assertTrue(ReportingDb.requireSummary(sampleId1, ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn(sampleId1), study1, coreCohort1));

        String sampleId2 = "CPCT01000001T";
        LimsCoreCohort coreCohort2 = LimsCoreCohort.fromSampleId(sampleId2);
        LimsStudy study2 = LimsStudy.fromSampleId(sampleId2);
        assertFalse(ReportingDb.requireSummary(sampleId2, ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn(sampleId2), study2, coreCohort2));

        String sampleId3 = "CORE01000001T";
        LimsCoreCohort coreCohort3 = LimsCoreCohort.fromSampleId(sampleId3);
        LimsStudy study3 = LimsStudy.fromSampleId(sampleId3);
        assertTrue(ReportingDb.requireSummary(sampleId3, ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn(sampleId3), study3, coreCohort3));

        String sampleId4 = "CORELR020000T";
        LimsCoreCohort coreCohort4 = LimsCoreCohort.fromSampleId(sampleId4);
        LimsStudy study4 = LimsStudy.fromSampleId(sampleId4);
        assertFalse(ReportingDb.requireSummary(sampleId4, ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn(sampleId4), study4, coreCohort4));

        String sampleId5 = "CORELR110000T";
        LimsCoreCohort coreCohort5 = LimsCoreCohort.fromSampleId(sampleId5);
        LimsStudy study5 = LimsStudy.fromSampleId(sampleId5);
        assertTrue(ReportingDb.requireSummary(sampleId5, ExampleAnalysisTestFactory.buildAnalysisWithAllTablesFilledIn(sampleId5), study5, coreCohort5));
    }
}