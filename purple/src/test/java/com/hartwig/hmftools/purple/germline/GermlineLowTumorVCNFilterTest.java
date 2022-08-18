package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.variant.VariantVcfTags.PURPLE_VARIANT_CN_INFO;
import static com.hartwig.hmftools.common.variant.Hotspot.HOTSPOT_FLAG;
import static com.hartwig.hmftools.common.variant.Hotspot.NEAR_HOTSPOT_FLAG;
import static com.hartwig.hmftools.purple.germline.GermlineLowTumorVCNFilter.MIN_QUAL_HOTSPOT;
import static com.hartwig.hmftools.purple.germline.GermlineLowTumorVCNFilter.MIN_QUAL_OTHER;
import static com.hartwig.hmftools.purple.germline.GermlineLowTumorVCNFilter.MIN_TUMOR_VCN;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.test.VariantContextFromString;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class GermlineLowTumorVCNFilterTest
{
    private static final double LOW_VN = MIN_TUMOR_VCN - 0.1;

    @Test
    public void testHighVariantCopyNumber()
    {
        assertFiltered(false, createVariant(MIN_QUAL_OTHER - 1, MIN_TUMOR_VCN, false));
    }

    @Test
    public void testHotspot()
    {
        assertFiltered(false, createVariant(MIN_QUAL_HOTSPOT, LOW_VN, true));
        assertFiltered(true, createVariant(MIN_QUAL_HOTSPOT - 1, LOW_VN, true));
    }

    @Test
    public void testOther()
    {
        assertFiltered(false, createVariant(MIN_QUAL_OTHER, LOW_VN, false));
        assertFiltered(true, createVariant(MIN_QUAL_OTHER - 1, LOW_VN, false));
        assertFiltered(true, createVariant(MIN_QUAL_HOTSPOT, LOW_VN, false));
    }

    private void assertFiltered(boolean expectedFiltered, VariantContext variantContext)
    {
        GermlineVariant variant = new GermlineVariant(variantContext);
        GermlineLowTumorVCNFilter.processVariant(variant);

        if(expectedFiltered)
        {
            assertTrue(variant.filters().contains(GermlineLowTumorVCNFilter.LOW_TUMOR_VCN_FILTER));
        }
        else
        {
            assertTrue(variant.isPass());
        }
    }

    @NotNull
    private static VariantContext createVariant(double qual, double variantCopyNumber, boolean isHotspot)
    {
        final String hotspotFlag = isHotspot ? HOTSPOT_FLAG : NEAR_HOTSPOT_FLAG;
        final String line = "11\t1000\tCOSM123;COSM456\tG\tA\t" + qual + "\tPASS\t" + PURPLE_VARIANT_CN_INFO + "=" + variantCopyNumber + ";"
                + hotspotFlag + "\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }
}
