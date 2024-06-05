package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.junit.Test;

public class GermlineVariantSelectorTest
{
    @Test
    public void canHandleNullGermlineVariants()
    {
        assertNull(GermlineVariantSelector.selectInterestingUnreportedVariants(null));
    }

    @Test
    public void canSelectInterestingUnreportedVariants()
    {
        PurpleVariant reportedHotspot = TestPurpleVariantFactory.builder().hotspot(HotspotType.HOTSPOT).build()
                .withCanonicalImpact(TestPurpleVariantFactory.impactBuilder().reported(true).build());
        PurpleVariant unreportedHotspot = TestPurpleVariantFactory.builder().hotspot(HotspotType.HOTSPOT).build();
        PurpleVariant unreportedNearHotspot = TestPurpleVariantFactory.builder().hotspot(HotspotType.NEAR_HOTSPOT).build();

        List<PurpleVariant> variants = Lists.newArrayList(reportedHotspot, unreportedHotspot, unreportedNearHotspot);
        List<PurpleVariant> interesting = GermlineVariantSelector.selectInterestingUnreportedVariants(variants);

        assertEquals(1, interesting.size());
        assertTrue(interesting.contains(unreportedHotspot));
    }
}