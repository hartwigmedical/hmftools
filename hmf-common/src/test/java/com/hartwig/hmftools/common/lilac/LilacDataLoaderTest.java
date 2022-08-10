package com.hartwig.hmftools.common.lilac;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.hla.LilacSummaryData;

import org.junit.Test;

public class LilacDataLoaderTest {

    private static final String LILAC_QC_CSV = Resources.getResource("lilac/sample.lilac.qc.csv").getPath();
    private static final String LILAC_RESULT_CSV = Resources.getResource("lilac/sample.lilac.csv").getPath();

    @Test
    public void canLoadFromTestFiles() throws IOException
    {
        LilacSummaryData summaryData = LilacSummaryData.load(LILAC_QC_CSV, LILAC_RESULT_CSV);
        assertNotNull(summaryData);
    }
}