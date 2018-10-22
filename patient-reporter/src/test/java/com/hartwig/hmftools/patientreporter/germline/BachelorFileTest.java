package com.hartwig.hmftools.patientreporter.germline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class BachelorFileTest {

    private static final double EPSILON = 1.0e-10;
    private static final String BACHELOR_FILE = Resources.getResource("germline/germline_variants.csv").getPath();

    @Test
    public void canReadTestBachelorFile() throws IOException {
        List<GermlineVariant> germlineVariants = BachelorFile.loadBachelorFile(BACHELOR_FILE);

        assertEquals(1, germlineVariants.size());

        GermlineVariant variant = germlineVariants.get(0);
        assertTrue(variant.passFilter());
        assertEquals("BRCA1", variant.gene());
        assertEquals("c.191G>A", variant.hgvsCodingImpact());
        assertEquals("p.Cys64Tyr", variant.hgvsProteinImpact());
        assertEquals(Strings.EMPTY, variant.germlineStatus());
        assertEquals(33, variant.alleleReadCount());
        assertEquals(89, variant.totalReadCount());
        assertEquals(2.99, variant.adjustedCopyNumber(), EPSILON);
        assertEquals(0.12, variant.adjustedVAF(), EPSILON);
        assertEquals(1.07, variant.minorAllelePloidy(), EPSILON);
        assertFalse(variant.biallelic());
    }
}