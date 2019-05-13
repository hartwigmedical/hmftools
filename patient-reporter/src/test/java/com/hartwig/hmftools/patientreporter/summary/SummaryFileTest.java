package com.hartwig.hmftools.patientreporter.summary;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class SummaryFileTest {

    private static final String SUMMARY_SAMPLES_CSV = Resources.getResource("csv/summary_samples.csv").getPath();

    @Test
    public void canLoadSummarySamplesCSV() throws IOException {
        SummaryModel summaryModel = SummaryFile.buildFromCsv(SUMMARY_SAMPLES_CSV);

        assertEquals(1, summaryModel.sizeSummarySamples().size());
    }
}