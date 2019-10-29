package com.hartwig.hmftools.patientreporter.reportingdb;

import static org.junit.Assert.assertEquals;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.ImmutableSampleMetadata;
import com.hartwig.hmftools.patientreporter.SampleMetadata;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReason;

import org.junit.Ignore;
import org.junit.Test;

public class ReportingDbTest {

    private static final String REPORT_DATES_TSV = Resources.getResource("lims/report_dates_lims.tsv").getPath();
    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final boolean WRITE_TO_TSV = false;

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity";
    private static final String PURPLE_QC = BASE_DIRECTORY + "/purple/sample.purple.qc";

    @Ignore
    @Test
    public void canReadReportDatesTsv() throws IOException {
        List<ReportingEntry> reportDates = ReportingDb.read(REPORT_DATES_TSV);

        assertEquals(2, reportDates.size());

        ReportingEntry reportingEntry1 = reportDates.get(0);
        assertEquals("sampleId", reportingEntry1.sampleId());
        assertEquals("ABCD", reportingEntry1.tumorBarcode());
        assertEquals("22/10/2019", reportingEntry1.reportDate());
        assertEquals("sequence_report", reportingEntry1.sourceReport());
        assertEquals("70%", reportingEntry1.purity());
        assertEquals("NORMAL", reportingEntry1.status());
        assertEquals("PASS", reportingEntry1.qcStatus());

        ReportingEntry reportingEntry2 = reportDates.get(1);
        assertEquals("sampleId", reportingEntry2.sampleId());
        assertEquals("EFGH", reportingEntry2.tumorBarcode());
        assertEquals("22/10/2019", reportingEntry2.reportDate());
        assertEquals("SHALLOW_SEQ_LOW_PURITY", reportingEntry2.sourceReport());
    }

    @Test
    public void canWriteReportDatesToTSV() throws IOException {
        if (WRITE_TO_TSV) {
            File reportDatesTsv = new File(REPORT_BASE_DIR + "/report_dates_lims.tsv");

            if (reportDatesTsv.createNewFile()) {
                BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTsv, true));
                writer.write("sampleId\ttumorBarcode\treportDate\tsourceReport\tpurity\tstatus\tqcStatus\n");
                writer.close();
            }

            SampleMetadata sampleMetaData = ImmutableSampleMetadata.builder()
                    .refSampleId("refSampleId")
                    .refSampleBarcode("EFGH")
                    .tumorSampleId("sampleId1")
                    .tumorSampleBarcode("ABCD")
                    .build();

            String clinicalSummary = "";

            ReportingDb.generateOutputReportDatesSeqRapports(reportDatesTsv.getPath(),
                    PURPLE_PURITY_TSV,
                    sampleMetaData,
                    PURPLE_QC,
                    false,
                    clinicalSummary,
                    true);

            QCFailReason reason = QCFailReason.fromIdentifier("shallow_seq_low_purity");

            ReportingDb.generateOutputReportDatesQCFailReport(reportDatesTsv.getPath(), reason, sampleMetaData, true);
        }
    }
}