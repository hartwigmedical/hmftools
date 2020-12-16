package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.patientreporter.PatientReportUtils;

import org.junit.Test;

public class SummaryFileTest {

    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();

    @Test
    public void summaryFromCSVWithNewLines() throws IOException {
        SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);
        assertEquals(1, summaryModel.summaryCount());

        LimsCohortConfig cohortConfig =
                PatientReportUtils.createCohortConfig("CORE", true, true, false, true, true, true, true, false, true, true, true, true);

        String summary = summaryModel.findSummaryForSample("sample", cohortConfig);

        assertEquals(3, summary.split("\n").length);
    }
}