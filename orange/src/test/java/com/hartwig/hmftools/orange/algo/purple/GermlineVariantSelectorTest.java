package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.Hotspot;

import org.junit.Test;

public class GermlineVariantSelectorTest {

    @Test
    public void canHandleNullGermlineVariants() {
        assertNull(GermlineVariantSelector.selectInterestingUnreportedVariants(null));
    }

    @Test
    public void canSelectInterestingUnreportedVariants() {
        PurpleVariant reportedHotspot = TestPurpleVariantFactory.builder().reported(true).hotspot(Hotspot.HOTSPOT).build();
        PurpleVariant unreportedHotspot = TestPurpleVariantFactory.builder().reported(false).hotspot(Hotspot.HOTSPOT).build();
        PurpleVariant unreportedNearHotspot = TestPurpleVariantFactory.builder().reported(false).hotspot(Hotspot.NEAR_HOTSPOT).build();

        List<PurpleVariant> variants = Lists.newArrayList(reportedHotspot, unreportedHotspot, unreportedNearHotspot);
        List<PurpleVariant> interesting = GermlineVariantSelector.selectInterestingUnreportedVariants(variants);

        assertEquals(1, interesting.size());
        assertTrue(interesting.contains(unreportedHotspot));
    }
}