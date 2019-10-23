package com.hartwig.hmftools.patientreporter.ReportDates;

import static org.junit.Assert.*;

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

public class ReportDatesAnalyzerTest {
    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();

    private static final String REPORT_DATES_TSV = Resources.getResource("lims/report_dates_lims.tsv").getPath();
    private static final String REPORT_BASE_DIR = System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp";
    private static final boolean WRITE_TO_TSV = true;

    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity";
    private static final String PURPLE_QC = BASE_DIRECTORY + "/purple/sample.purple.qc";

    @Ignore
    @Test
    public void canReadReportDatesTsv() throws IOException {
        List<ReportDates> reportDates = ReportDatesAnalyzer.read(REPORT_DATES_TSV);

        assertEquals(2, reportDates.size());

        ReportDates reportDates1 = reportDates.get(0);
        assertEquals("sampleId", reportDates1.sampleId());
        assertEquals("ABCD", reportDates1.tumorBarcode());
        assertEquals("22/10/2019", reportDates1.reportDate());
        assertEquals("sequence_report", reportDates1.sourceReport());
        assertEquals("70%", reportDates1.purity());
        assertEquals("NORMAL", reportDates1.status());
        assertEquals("PASS", reportDates1.qcStatus());

        ReportDates reportDates2 = reportDates.get(1);
        assertEquals("sampleId", reportDates2.sampleId());
        assertEquals("EFGH", reportDates2.tumorBarcode());
        assertEquals("22/10/2019", reportDates2.reportDate());
        assertEquals("SHALLOW_SEQ_LOW_PURITY", reportDates2.sourceReport());
    }

    @Test
    public void canWriteReportDatesToTSV() throws IOException {
        if (WRITE_TO_TSV) {
            File reportDatesTSV = new File(REPORT_BASE_DIR + "/report_dates_lims.tsv");
            reportDatesTSV.createNewFile();

            if (reportDatesTSV.length() == 0) {
                BufferedWriter writer = new BufferedWriter(new FileWriter(reportDatesTSV, true));
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

            ReportDatesAnalyzer.generateOutputReportDatesSeqRapports(reportDatesTSV.getPath(),
                    PURPLE_PURITY_TSV,
                    sampleMetaData,
                    PURPLE_QC,
                    false,
                    clinicalSummary,
                    WRITE_TO_TSV);

            QCFailReason reason = QCFailReason.fromIdentifier("shallow_seq_low_purity");

            ReportDatesAnalyzer.generateOutputReportDatesQCFailReport(reason, reportDatesTSV.getPath(), sampleMetaData, WRITE_TO_TSV);
        }
    }
}