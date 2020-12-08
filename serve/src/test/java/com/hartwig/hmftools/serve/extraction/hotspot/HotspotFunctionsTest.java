package com.hartwig.hmftools.serve.extraction.hotspot;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HotspotFunctionsTest {

    @Test
    public void canConsolidateEmptyHotspots() {
        assertTrue(HotspotFunctions.consolidate(Lists.newArrayList()).isEmpty());
    }

    @Test
    public void canConsolidateHotspotsFromOneSource() {
        Knowledgebase source = Knowledgebase.HARTWIG_CURATED;
        Set<KnownHotspot> knownHotspots = Sets.newHashSet();
        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot1())
                .addSources(source)
                .gene("gene1")
                .transcript("trans1")
                .proteinAnnotation("prot1")
                .build());

        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot1())
                .addSources(source)
                .gene("gene1")
                .transcript(null)
                .proteinAnnotation("prot2")
                .build());

        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot2())
                .addSources(source)
                .gene("gene2")
                .transcript("trans2")
                .proteinAnnotation("prot3")
                .build());

        Set<KnownHotspot> consolidateHotspots = HotspotFunctions.consolidate(knownHotspots);
        assertEquals(2, consolidateHotspots.size());

        KnownHotspot gene1 = findByGene(consolidateHotspots, "gene1");
        assertEquals(Sets.newHashSet(source), gene1.sources());
        assertEquals("trans1", gene1.transcript());
        assertEquals("prot1", gene1.proteinAnnotation());

        KnownHotspot gene2 = findByGene(consolidateHotspots, "gene2");
        assertEquals(Sets.newHashSet(source), gene2.sources());
        assertEquals("trans2", gene2.transcript());
        assertEquals("prot3", gene2.proteinAnnotation());
    }

    @Test
    public void canConsolidateHotspotsFromTwoSources() {
        String gene = "gene1";
        Knowledgebase source1 = Knowledgebase.HARTWIG_CURATED;
        Knowledgebase source2 = Knowledgebase.HARTWIG_COHORT;
        Set<KnownHotspot> knownHotspots = Sets.newHashSet();
        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot1())
                .addSources(source1)
                .gene(gene)
                .transcript(null)
                .proteinAnnotation("prot1")
                .build());

        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot1())
                .addSources(source2)
                .gene(gene)
                .transcript("trans2")
                .proteinAnnotation("prot2")
                .build());

        Set<KnownHotspot> consolidateHotspots = HotspotFunctions.consolidate(knownHotspots);
        assertEquals(1, consolidateHotspots.size());

        KnownHotspot hotspot = findByGene(consolidateHotspots, gene);
        assertEquals(Sets.newHashSet(source1, source2), hotspot.sources());
        assertEquals("trans2", hotspot.transcript());
        assertEquals("prot2", hotspot.proteinAnnotation());
    }

    @NotNull
    private static VariantHotspot hotspot1() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build();
    }

    @NotNull
    private static VariantHotspot hotspot2() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").position(20).ref("A").alt("T").build();
    }

    @NotNull
    private static KnownHotspot findByGene(@NotNull Iterable<KnownHotspot> hotspots, @NotNull String gene) {
        for (KnownHotspot hotspot : hotspots) {
            if (hotspot.gene().equals(gene)) {
                return hotspot;
            }
        }

        throw new IllegalStateException("Could not find gene in hotspots: " + gene);
    }
}