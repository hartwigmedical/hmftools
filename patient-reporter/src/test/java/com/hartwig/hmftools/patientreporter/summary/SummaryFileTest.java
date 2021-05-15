package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortConfig;
import com.hartwig.hmftools.common.lims.cohort.LimsCohortTestFactory;

import org.junit.Test;

public class SummaryFileTest {

    private static final String SAMPLE_SUMMARY_TSV = Resources.getResource("sample_summary/sample_summary.tsv").getPath();

    @Test
    public void summaryFromTsvWithNewLines() throws IOException {
        SummaryModel summaryModel = SummaryFile.buildFromTsv(SAMPLE_SUMMARY_TSV);
        assertEquals(1, summaryModel.summaryCount());

        LimsCohortConfig cohortConfig = LimsCohortTestFactory.createCOLOCohortConfig();

        String summary = summaryModel.findSummaryForSample("sample", cohortConfig);

        assertEquals(3, summary.split("\n").length);
    }
}