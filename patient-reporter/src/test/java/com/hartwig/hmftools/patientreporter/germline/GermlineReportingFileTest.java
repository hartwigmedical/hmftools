package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class GermlineReportingFileTest {

    private static final String GERMLINE_REPORTING_TSV = Resources.getResource("germline_reporting/germline_reporting.tsv").getPath();

    @Test
    public void canLoadGermlineReportingTsv() throws IOException {
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(GERMLINE_REPORTING_TSV);

        assertEquals(4, germlineReportingModel.entries().size());

        GermlineReportingEntry apc = germlineReportingModel.entryForGene("APC");
        assertNotNull(apc);
        assertEquals(GermlineCondition.ALWAYS, apc.condition());
        assertNull(apc.conditionFilter());
    }
}