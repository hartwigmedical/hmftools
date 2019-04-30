package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.*;

import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineGenesReporting;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineGenesReportingFile;

import org.junit.Test;

public class GermlineGenesReportingFileTest {

    private static final String GERMLINE_GENES_REPORTING_CSV = Resources.getResource("csv/germline_genes_reporting.csv").getPath();

    @Test
    public void canLoadGermlineGenesReportingCsv() throws IOException {
        GermlineGenesReporting germlineGenesReporting = GermlineGenesReportingFile.buildFromCsv(GERMLINE_GENES_REPORTING_CSV);

        assertEquals(3, germlineGenesReporting.germlineGenes().size());
    }
}