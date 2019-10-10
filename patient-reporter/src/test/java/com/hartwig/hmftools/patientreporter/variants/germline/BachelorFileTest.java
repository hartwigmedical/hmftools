package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class BachelorFileTest {

    private static final double EPSILON = 1.0e-10;
    private static final String BACHELOR_FILE = Resources.getResource("test_run/bachelor/sample_germline_variants.csv").getPath();

    @Test
    public void canReadTestBachelorFile() throws IOException {
        List<GermlineVariant> germlineVariants = BachelorFile.loadBachelorCsv(BACHELOR_FILE);

        assertEquals(1, germlineVariants.size());

        GermlineVariant variant = germlineVariants.get(0);
        assertTrue(variant.passFilter());
        assertEquals("BRCA1", variant.gene());
        assertEquals("c.68_69delAG", variant.hgvsCodingImpact());
        assertEquals("p.Glu23fs", variant.hgvsProteinImpact());
        assertEquals(45, variant.alleleReadCount());
        assertEquals(97, variant.totalReadCount());
        assertEquals(3.8773, variant.adjustedCopyNumber(), EPSILON);
        assertEquals(0.43355030873174055, variant.adjustedVAF(), EPSILON);
        assertFalse(variant.biallelic());
    }
}