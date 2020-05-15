package com.hartwig.hmftools.patientreporter.variants.germline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.junit.Test;

public class BachelorFileTest {

    private static final double EPSILON = 1.0e-10;
    private static final String BACHELOR_TSV = Resources.getResource("test_run/bachelor/sample.bachelor.germline_variant.tsv").getPath();

    @Test
    public void canReadTestBachelorFile() throws IOException {
        List<GermlineVariant> germlineVariants = BachelorFile.loadBachelorTsv(BACHELOR_TSV);

        assertEquals(1, germlineVariants.size());

        GermlineVariant variant = germlineVariants.get(0);
        assertEquals("17", variant.chromosome());
        assertEquals(41276044, variant.position());
        assertTrue(variant.passFilter());
        assertEquals("ACT", variant.ref());
        assertEquals("A", variant.alt());
        assertEquals("BRCA1", variant.gene());
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, variant.codingEffect());
        assertEquals("c.68_69delAG", variant.hgvsCodingImpact());
        assertEquals("p.Glu23fs", variant.hgvsProteinImpact());
        assertEquals(45, variant.alleleReadCount());
        assertEquals(97, variant.totalReadCount());
        assertEquals(0.43355030873174055, variant.adjustedVAF(), EPSILON);
        assertEquals(3.8773, variant.adjustedCopyNumber(), EPSILON);
        assertFalse(variant.biallelic());
    }
}