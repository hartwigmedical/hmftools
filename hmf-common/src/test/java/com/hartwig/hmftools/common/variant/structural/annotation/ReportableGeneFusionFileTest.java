package com.hartwig.hmftools.common.variant.structural.annotation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class ReportableGeneFusionFileTest {

    /*
    private static final String DETAILED_FUSION_FILE =Resources.getResource("variant/structural/fusions_detailed.csv").getPath();

    @Test
    public void canReadDetailedFusionFile() throws IOException
    {
        List<ReportableGeneFusion> fusions = ReportableGeneFusionFile.readFromDetailedFusionFile(DETAILED_FUSION_FILE);

        assertEquals(1, fusions.size());

        ReportableGeneFusion fusion = fusions.get(0);

        assertEquals("GEMIN5", fusion.geneStart());
        assertEquals("BRAF", fusion.geneEnd());
        assertEquals("ENST00000285873", fusion.geneTranscriptStart());
        assertEquals("ENST00000288602", fusion.geneTranscriptEnd());
        assertEquals("Intron 24", fusion.geneContextStart());
        assertEquals("Intron 8", fusion.geneContextEnd());
        Double ploidy = fusion.ploidy();
        assertNotNull(ploidy);
        assertEquals(1.98, ploidy, 1.0E-10);
    }
    */
}