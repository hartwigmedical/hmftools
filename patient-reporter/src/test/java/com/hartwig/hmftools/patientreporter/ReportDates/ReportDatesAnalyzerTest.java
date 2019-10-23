package com.hartwig.hmftools.patientreporter.ReportDates;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class ReportDatesAnalyzerTest {
    private static final String REPORT_DATES_TSV = Resources.getResource("lims/report_dates_lims.tsv").getPath();

    @Test
    public void canReadReportDatesTsv() throws IOException {
        List<ReportDates> reportDates = ReportDatesAnalyzer.read(REPORT_DATES_TSV);

        assertEquals(3, reportDates.size());

        ReportDates reportDates1 = reportDates.get(1);
        assertEquals("sampleId", reportDates1.sampleId());
        assertEquals("ABCD", reportDates1.tumorBarcode());
        assertEquals("22/10/2019", reportDates1.reportDate());
        assertEquals("sequence_report", reportDates1.sourceReport());
        assertEquals("70%", reportDates1.purity());
        assertEquals("NORMAL", reportDates1.status());
        assertEquals("PASS", reportDates1.qcStatus());

        ReportDates reportDates2 = reportDates.get(2);
        assertEquals("sampleId", reportDates2.sampleId());
        assertEquals("EFGH", reportDates2.tumorBarcode());
        assertEquals("22/10/2019", reportDates2.reportDate());
        assertEquals("SHALLOW_SEQ_LOW_PURITY", reportDates2.sourceReport());


    }

}