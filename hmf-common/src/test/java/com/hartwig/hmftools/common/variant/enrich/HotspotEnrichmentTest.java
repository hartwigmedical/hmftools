package com.hartwig.hmftools.common.variant.enrich;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.VariantContextFromString;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class HotspotEnrichmentTest {

    @Test
    public void testOverlap() {
        final String ponRef = "GATTACA";
        final String variantRef = "T";

        final VariantHotspot hotspot = ImmutableVariantHotspot.builder().chromosome("11").position(100).ref(ponRef).alt("A").build();
        // Max should be 100 + 7 - 1 + 5
        // Min should be 100 - 1 + 1 - 5

        assertOverlap(false, hotspot, 93, variantRef);
        assertOverlap(false, hotspot, 94, variantRef);
        assertOverlap(true, hotspot, 95, variantRef);
        assertOverlap(true, hotspot, 96, variantRef);

        assertOverlap(true, hotspot, 99, variantRef);
        assertOverlap(true, hotspot, 100, variantRef);
        assertOverlap(true, hotspot, 101, variantRef);

        assertOverlap(true, hotspot, 110, variantRef);
        assertOverlap(true, hotspot, 111, variantRef);
        assertOverlap(false, hotspot, 112, variantRef);
        assertOverlap(false, hotspot, 113, variantRef);
    }

    private void assertOverlap(boolean expected, @NotNull VariantHotspot hotspot, int variantStart, @NotNull final String variantRef) {
        final VariantContext variant = createNonHotspot(variantStart, variantRef);
        assertEquals(expected, HotspotEnrichment.overlaps(hotspot, variant));
    }

    @NotNull
    final VariantContext createNonHotspot(int start, @NotNull final String ref) {

        final String line = "11\t" + start + "\tCOSM123;COSM456\t" + ref
                + "\tC\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);

    }

}
