package com.hartwig.hmftools.purple.somatic;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.HotspotType;
import com.hartwig.hmftools.common.test.VariantContextFromString;
import com.hartwig.hmftools.common.variant.SimpleVariant;

import org.jetbrains.annotations.NotNull;
import org.junit.Assert;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class HotspotEnrichmentTest
{
    @Test
    public void testOverlap()
    {
        final String hotspotRef = "GATTACA";
        final String variantRef = "T";

        final SimpleVariant hotspot = new SimpleVariant("11", 100, hotspotRef, "A");
        // Max should be 100 + 7 - 1 + 5
        // Min should be 100 - 1 + 1 - 5

        assertOverlap(HotspotType.NON_HOTSPOT, hotspot, 93, variantRef);
        assertOverlap(HotspotType.NON_HOTSPOT, hotspot, 94, variantRef);
        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 95, variantRef);
        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 96, variantRef);

        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 99, variantRef);
        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 100, variantRef);
        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 101, variantRef);

        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 110, variantRef);
        assertOverlap(HotspotType.NEAR_HOTSPOT, hotspot, 111, variantRef);
        assertOverlap(HotspotType.NON_HOTSPOT, hotspot, 112, variantRef);
        assertOverlap(HotspotType.NON_HOTSPOT, hotspot, 113, variantRef);
    }

    @Test
    public void testExactMatch()
    {
        String hotspotRef = "GATTACA";

        SimpleVariant hotspot = new SimpleVariant("11", 100, hotspotRef, "A");
        assertOverlap(HotspotType.HOTSPOT, hotspot, 100, hotspotRef);
    }

    @Test
    public void testFromVariant()
    {
        VariantContext variant = createNonHotspotV37(1, "G");
        VariantContext nonHotspot = new VariantContextBuilder(variant).attribute("HOTSPOT", false).make();
        Assert.assertEquals(HotspotType.NON_HOTSPOT, HotspotType.fromVariant(nonHotspot));
    }

    private static void assertOverlap(HotspotType expected, SimpleVariant hotspot, int variantStart, String variantRef)
    {
        ListMultimap<Chromosome,SimpleVariant> hotspotMap = ArrayListMultimap.create();
        hotspotMap.put(HumanChromosome.fromString(hotspot.chromosome()), hotspot);
        HotspotEnrichment enrichment = new HotspotEnrichment(hotspotMap, true);

        VariantContext v37Variant = createNonHotspotV37(variantStart, variantRef);
        VariantContext v38Variant = createNonHotspotV38(variantStart, variantRef);

        enrichment.processVariant(v37Variant);
        enrichment.processVariant(v38Variant);

        assertEquals(expected, HotspotType.fromVariant(v37Variant));
        assertEquals(expected, HotspotType.fromVariant(v38Variant));
    }

    @NotNull
    private static VariantContext createNonHotspotV37(int start, final String ref)
    {
        final String line = "11\t" + start + "\tCOSM123;COSM456\t" + ref
                + "\tA\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }

    @NotNull
    private static VariantContext createNonHotspotV38(int start, final String ref)
    {
        final String line = "chr11\t" + start + "\tCOSM123;COSM456\t" + ref
                + "\tA\t.\tPASS\tCOSM2ENST=COSM123|GENE_TRANS1|c.1A>G|p.E1E|1,COSM456|GENE_TRANS2|c.2A>G|p.E2E|1\tGT:AD:DP\t0/1:73,17:91";
        return VariantContextFromString.decode(line);
    }
}
