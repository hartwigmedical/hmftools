package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.purple.germline.GermlineGenotypeEnrichment.GermlineGenotypeStatus.HET;
import static com.hartwig.hmftools.purple.germline.GermlineGenotypeEnrichment.GermlineGenotypeStatus.HOM_ALT;
import static com.hartwig.hmftools.purple.germline.GermlineGenotypeEnrichment.GermlineGenotypeStatus.LOW_VAF;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import org.junit.Test;

public class GermlineGenotypeEnrichmentTest
{
    @Test
    public void testCombinedStatus()
    {
        assertEquals(LOW_VAF, GermlineGenotypeEnrichment.combined(LOW_VAF, LOW_VAF));
        assertEquals(LOW_VAF, GermlineGenotypeEnrichment.combined(LOW_VAF, HET));
        assertEquals(LOW_VAF, GermlineGenotypeEnrichment.combined(LOW_VAF, HOM_ALT));

        assertEquals(HET, GermlineGenotypeEnrichment.combined(HET, LOW_VAF));
        assertEquals(HET, GermlineGenotypeEnrichment.combined(HET, HET));
        assertEquals(HET, GermlineGenotypeEnrichment.combined(HET, HOM_ALT));

        assertEquals(HET, GermlineGenotypeEnrichment.combined(HOM_ALT, LOW_VAF));
        assertEquals(HET, GermlineGenotypeEnrichment.combined(HOM_ALT, HET));
        assertEquals(HOM_ALT, GermlineGenotypeEnrichment.combined(HOM_ALT, HOM_ALT));
    }

    @Test
    public void testStatus()
    {
        assertStatus(LOW_VAF, 1, 17);
        assertStatus(LOW_VAF, 2, 21);

        assertStatus(HET, 1, 16);
        assertStatus(HET, 2, 20);

        assertStatus(HET, 13, 14);
        assertStatus(HOM_ALT, 14, 15);
    }

    private static void assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus expected, int alleleReadCount, int totalReadCount)
    {
        AllelicDepth depth = new AllelicDepth(totalReadCount, alleleReadCount);
        assertEquals(expected, GermlineGenotypeEnrichment.status(depth));
    }
}
