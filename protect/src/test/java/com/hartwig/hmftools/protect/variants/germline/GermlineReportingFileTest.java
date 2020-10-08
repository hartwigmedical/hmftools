package com.hartwig.hmftools.protect.variants.germline;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class GermlineReportingFileTest {

    private static final String GERMLINE_GENES_REPORTING_CSV = Resources.getResource("germline/germline_genes_reporting.csv").getPath();

    @Test
    public void canLoadGermlineGenesReportingCsv() throws IOException {
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(GERMLINE_GENES_REPORTING_CSV);

        assertEquals(3, germlineReportingModel.reportableGermlineGenes().size());
    }
}