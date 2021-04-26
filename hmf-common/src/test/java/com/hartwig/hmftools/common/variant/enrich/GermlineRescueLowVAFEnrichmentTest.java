package com.hartwig.hmftools.common.variant.enrich;

import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_VARIANT_CN_INFO;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.VariantContextFromString;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineRescueLowVAFEnrichmentTest {

    @Test
    public void testRestore() {
        assertFiltered(true, createVariant(0.4, 2));
        assertFiltered(true, createVariant(0.5, 2));
        assertFiltered(false, createVariant(0.5, 3));
    }

    private static void assertFiltered(boolean expectedFiltered, VariantContext victim) {
        VariantContext updated = GermlineRescueLowVAFEnrichment.process(VariantContextFromString.SAMPLE, victim);
        assertEquals(expectedFiltered, updated.isFiltered());
        if (expectedFiltered) {
            assertTrue(updated.getFilters().contains(GermlineGenotypeEnrichment.LOW_VAF_FILTER));
        }
    }

    @NotNull
    private static VariantContext createVariant(double variantCopyNumber, int alleleReadCount) {
        final String line =
                "11\t1000\tCOSM123;COSM456\tG\tA\t100\t" + GermlineGenotypeEnrichment.LOW_VAF_FILTER + "\t" + PURPLE_VARIANT_CN_INFO + "="
                        + variantCopyNumber + ";NEAR_HOTSPOT\tGT:AD:DP\t0/1:73," + alleleReadCount + ":91";
        return VariantContextFromString.decode(line);
    }
}
