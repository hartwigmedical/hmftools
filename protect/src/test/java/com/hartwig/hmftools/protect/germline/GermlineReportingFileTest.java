package com.hartwig.hmftools.protect.germline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class GermlineReportingFileTest {

    private static final String GERMLINE_REPORTING_TSV = Resources.getResource("germline/germline_reporting.tsv").getPath();

    @Test
    public void canLoadGermlineReportingTsv() throws IOException {
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(GERMLINE_REPORTING_TSV);

        assertEquals(3, germlineReportingModel.entries().size());

        GermlineReportingEntry apc = germlineReportingModel.entryForGene("APC");
        assertNotNull(apc);
        assertTrue(apc.notifyClinicalGeneticist());
        assertNull(apc.exclusiveHgvsProteinFilter());
    }
}