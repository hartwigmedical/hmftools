package com.hartwig.hmftools.serve.hotspot;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

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
        List<KnownHotspot> knownHotspots = Lists.newArrayList();
        assertTrue(HotspotFunctions.consolidate(knownHotspots).isEmpty());
    }

    @Test
    public void canConsolidateHotspotsFromOneSource() {
        Knowledgebase source = Knowledgebase.HARTWIG_CURATED;
        List<KnownHotspot> knownHotspots = Lists.newArrayList();
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

        List<KnownHotspot> consolidateHotspots = HotspotFunctions.consolidate(knownHotspots);
        assertEquals(2, consolidateHotspots.size());

        assertEquals(Sets.newHashSet(source), consolidateHotspots.get(0).sources());
        assertEquals("gene1", consolidateHotspots.get(0).gene());
        assertEquals("trans1", consolidateHotspots.get(0).transcript());
        assertEquals("prot1", consolidateHotspots.get(0).proteinAnnotation());

        assertEquals(Sets.newHashSet(source), consolidateHotspots.get(1).sources());
        assertEquals("gene2", consolidateHotspots.get(1).gene());
        assertEquals("trans2", consolidateHotspots.get(1).transcript());
        assertEquals("prot3", consolidateHotspots.get(1).proteinAnnotation());
    }

    @Test
    public void canConsolidateHotspotsFromTwoSources() {
        Knowledgebase source1 = Knowledgebase.HARTWIG_CURATED;
        Knowledgebase source2 = Knowledgebase.HARTWIG_COHORT;
        List<KnownHotspot> knownHotspots = Lists.newArrayList();
        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot1())
                .addSources(source1)
                .gene("gene1")
                .transcript(null)
                .proteinAnnotation("prot1")
                .build());

        knownHotspots.add(ImmutableKnownHotspot.builder()
                .from(hotspot1())
                .addSources(source2)
                .gene("gene1")
                .transcript("trans2")
                .proteinAnnotation("prot2")
                .build());

        List<KnownHotspot> consolidateHotspots = HotspotFunctions.consolidate(knownHotspots);
        assertEquals(1, consolidateHotspots.size());

        assertEquals(Sets.newHashSet(source1, source2), consolidateHotspots.get(0).sources());
        assertEquals("gene1", consolidateHotspots.get(0).gene());
        assertEquals("trans2", consolidateHotspots.get(0).transcript());
        assertEquals("prot2", consolidateHotspots.get(0).proteinAnnotation());
    }

    @NotNull
    private static VariantHotspot hotspot1() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build();
    }

    @NotNull
    private static VariantHotspot hotspot2() {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").position(20).ref("A").alt("T").build();
    }
}