package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.VariantVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_VARIANT_CN_INFO;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.VariantContextFromString;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineRescueLowVAFTest
{
    @Test
    public void testRestore()
    {
        assertFiltered(true, createVariant(0.4, 2));
        assertFiltered(true, createVariant(0.5, 2));
        assertFiltered(false, createVariant(0.5, 3));
    }

    private static void assertFiltered(boolean expectedFiltered, VariantContext variantContext)
    {
        GermlineRescueLowVAF germlineRescueLowVAF = new GermlineRescueLowVAF(VariantContextFromString.SAMPLE);

        GermlineVariant variant = new GermlineVariant(variantContext);
        germlineRescueLowVAF.processVariant(variant);

        if(expectedFiltered)
        {
            assertTrue(variant.filters().contains(GermlineGenotypeEnrichment.LOW_VAF_FILTER));
        }
        else
        {
            assertTrue(variant.isPass());
        }
    }

    @NotNull
    private static VariantContext createVariant(double variantCopyNumber, int alleleReadCount)
    {
        final String line =
                "11\t1000\tCOSM123;COSM456\tG\tA\t100\t" + GermlineGenotypeEnrichment.LOW_VAF_FILTER + "\t" + PURPLE_VARIANT_CN_INFO + "="
                        + variantCopyNumber + ";NEAR_HOTSPOT\tGT:AD:DP\t0/1:73," + alleleReadCount + ":91";

        return VariantContextFromString.decode(line);
    }
}
