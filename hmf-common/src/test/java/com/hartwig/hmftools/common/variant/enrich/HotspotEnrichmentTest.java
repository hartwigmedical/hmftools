package com.hartwig.hmftools.common.variant.enrich;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantContextFromString;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class HotspotEnrichmentTest {

    @Test
    public void testOverlap() {
        final String hotspotRef = "GATTACA";
        final String variantRef = "T";

        final VariantHotspot hotspot =
                ImmutableVariantHotspotImpl.builder().chromosome("11").position(100).ref(hotspotRef).alt("A").build();
        // Max should be 100 + 7 - 1 + 5
        // Min should be 100 - 1 + 1 - 5

        assertOverlap(Hotspot.NON_HOTSPOT, hotspot, 93, variantRef);
        assertOverlap(Hotspot.NON_HOTSPOT, hotspot, 94, variantRef);
        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 95, variantRef);
        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 96, variantRef);

        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 99, variantRef);
        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 100, variantRef);
        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 101, variantRef);

        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 110, variantRef);
        assertOverlap(Hotspot.NEAR_HOTSPOT, hotspot, 111, variantRef);
        assertOverlap(Hotspot.NON_HOTSPOT, hotspot, 112, variantRef);
        assertOverlap(Hotspot.NON_HOTSPOT, hotspot, 113, variantRef);
    }

    @Test
    public void testExactMatch() {
        final String hotspotRef = "GATTACA";

        final VariantHotspot hotspot =
                ImmutableVariantHotspotImpl.builder().chromosome("11").position(100).ref(hotspotRef).alt("A").build();
        assertOverlap(Hotspot.HOTSPOT, hotspot, 100, hotspotRef);
    }

    @Test
    public void testFromVariant() {
        VariantContext variant = createNonHotspotV37(1, "G");
        VariantContext nonHotspot = new VariantContextBuilder(variant).attribute("HOTSPOT", false).make();
        assertEquals(Hotspot.NON_HOTSPOT, HotspotEnrichment.fromVariant(nonHotspot));
    }

    private static void assertOverlap(@NotNull Hotspot expected, @NotNull VariantHotspot hotspot, int variantStart,
            @NotNull String variantRef) {
        final List<VariantContext> result = Lists.newArrayList();

        final VariantHotspotEnrichment enrichment =
                new VariantHotspotEnrichment(Multimaps.fromPositions(Lists.newArrayList(hotspot)), result::add);

        final VariantContext v37Variant = createNonHotspotV37(variantStart, variantRef);
        final VariantContext v38Variant = createNonHotspotV38(variantStart, variantRef);

        enrichment.accept(v37Variant);
        enrichment.accept(v38Variant);

        assertEquals(expected, HotspotEnrichment.fromVariant(result.get(0)));
        assertEquals(expected, HotspotEnrichment.fromVariant(result.get(1)));
    }

    @NotNull
    private static VariantContext createNonHotspotV37(int start, @NotNull final String ref) {
        final String line = "11\t" + start + "\tCOSM123;COSM456\t" + ref
                + "\tA\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }

    @NotNull
    private static VariantContext createNonHotspotV38(int start, @NotNull final String ref) {
        final String line = "chr11\t" + start + "\tCOSM123;COSM456\t" + ref
                + "\tA\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }
}
