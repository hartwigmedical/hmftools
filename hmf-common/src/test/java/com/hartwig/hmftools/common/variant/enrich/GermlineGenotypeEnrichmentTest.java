package com.hartwig.hmftools.common.variant.enrich;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineGenotypeEnrichmentTest {

    @Test
    public void testStatus() {
        assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus.LOW_VAF, 1, 17);
        assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus.LOW_VAF, 2, 21);

        assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus.HET, 1, 16);
        assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus.HET, 2, 20);

        assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus.HET, 13, 14);
        assertStatus(GermlineGenotypeEnrichment.GermlineGenotypeStatus.HOM, 14, 15);
    }

    private static void assertStatus(@NotNull GermlineGenotypeEnrichment.GermlineGenotypeStatus expected, int alleleReadCount,
            int totalReadCount) {
        AllelicDepth depth = ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
        assertEquals(expected, GermlineGenotypeEnrichment.status(depth));
    }
}
