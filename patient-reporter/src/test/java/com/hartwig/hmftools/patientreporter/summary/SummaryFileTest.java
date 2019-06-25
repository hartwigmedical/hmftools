package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class SummaryFileTest {

    private static final String SAMPLE_SUMMARY_CSV = Resources.getResource("csv/sample_summary.csv").getPath();

    @Test
    public void canLoadSummarySamplesCSV() throws IOException {
        SummaryModel summaryModel = SummaryFile.buildFromCsv(SAMPLE_SUMMARY_CSV);

        assertEquals(1, summaryModel.summaryCount());
    }
}