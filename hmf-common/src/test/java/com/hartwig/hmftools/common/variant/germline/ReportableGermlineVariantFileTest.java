package com.hartwig.hmftools.common.variant.germline;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.junit.Test;

public class ReportableGermlineVariantFileTest {

    private static final double EPSILON = 1.0e-10;
    private static final String BACHELOR_TSV = Resources.getResource("variant/germline/sample.reportable_germline_variant.tsv").getPath();

    @Test
    public void canReadTestBachelorFile() throws IOException {
        List<ReportableGermlineVariant> germlineVariants = ReportableGermlineVariantFile.read(BACHELOR_TSV);

        assertEquals(1, germlineVariants.size());

        ReportableGermlineVariant variant = germlineVariants.get(0);
        assertEquals("BRCA1", variant.gene());
        assertEquals("17", variant.chromosome());
        assertEquals(41276044, variant.position());
        assertEquals("ACT", variant.ref());
        assertEquals("A", variant.alt());
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, variant.codingEffect());
        assertEquals("c.68_69delAG", variant.hgvsCoding());
        assertEquals("p.Glu23fs", variant.hgvsProtein());
        assertEquals(45, variant.alleleReadCount());
        assertEquals(97, variant.totalReadCount());
        assertEquals(0.43355030873174055, variant.adjustedVaf(), EPSILON);
        assertEquals(3.8773, variant.adjustedCopyNumber(), EPSILON);
        assertFalse(variant.biallelic());
    }

}