package com.hartwig.hmftools.serve.hotspot;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class HotspotFunctionsTest {

    @Test
    public void canConvertHotspots() {
        TestHotspotSourceEntry entry1 = new TestHotspotSourceEntry("gene1", null, "annotation1");
        TestHotspotSourceEntry entry2 = new TestHotspotSourceEntry("gene1", "transcript1", "annotation2");
        Map<TestHotspotSourceEntry, List<VariantHotspot>> map = Maps.newHashMap();

        VariantHotspot hotspot1 = ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("AA").alt("TT").build();
        VariantHotspot hotspot2 = ImmutableVariantHotspotImpl.builder().chromosome("1").position(11).ref("A").alt("T").build();
        VariantHotspot hotspot3 = ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("AA").alt("TT").build();

        map.put(entry1, Lists.newArrayList(hotspot1, hotspot2));
        map.put(entry2, Lists.newArrayList(hotspot3));

        String source = "source";
        Map<VariantHotspot, HotspotAnnotation> converted = HotspotFunctions.convertHotspotMap(source, map);

        assertEquals(2, converted.size());
        HotspotAnnotation annotation1 = converted.get(hotspot1);
        assertEquals("gene1", annotation1.gene());
        assertEquals("transcript1", annotation1.transcript());
        assertEquals("annotation2", annotation1.proteinAnnotation());
        assertEquals(Sets.newHashSet(source), annotation1.sources());

        HotspotAnnotation annotation2 = converted.get(hotspot2);
        assertEquals("gene1", annotation2.gene());
        assertNull(annotation2.transcript());
        assertEquals("annotation1", annotation2.proteinAnnotation());
        assertEquals(Sets.newHashSet(source), annotation2.sources());
    }

    @Test
    public void canMergeHotspots() {
        String source1 = "source1";
        String source2 = "source2";

        VariantHotspot hotspot1 = ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build();
        VariantHotspot hotspot2 = ImmutableVariantHotspotImpl.builder().chromosome("1").position(20).ref("A").alt("T").build();
        HotspotAnnotation source1Annotation1 = new HotspotAnnotation(Sets.newHashSet(source1), "gene1", "transcript1", "annotation1");
        HotspotAnnotation source1Annotation2 = new HotspotAnnotation(Sets.newHashSet(source1), "gene1", "transcript1", "annotation2");
        HotspotAnnotation source2Annotation1 = new HotspotAnnotation(Sets.newHashSet(source2), "gene1", "transcript1", "annotation1");

        Map<VariantHotspot, HotspotAnnotation> source1Map = Maps.newHashMap();
        source1Map.put(hotspot1, source1Annotation1);
        source1Map.put(hotspot2, source1Annotation2);

        Map<VariantHotspot, HotspotAnnotation> source2Map = Maps.newHashMap();
        source2Map.put(hotspot1, source2Annotation1);

        Map<VariantHotspot, HotspotAnnotation> mergedMap = HotspotFunctions.mergeHotspotMaps(Lists.newArrayList(source1Map, source2Map));

        assertEquals(2, mergedMap.size());
        HotspotAnnotation mergedAnnotation1 = mergedMap.get(hotspot1);
        assertEquals(Sets.newHashSet(source1, source2), mergedAnnotation1.sources());

        HotspotAnnotation mergedAnnotation2 = mergedMap.get(hotspot2);
        assertEquals(Sets.newHashSet(source1), mergedAnnotation2.sources());
    }

    private static class TestHotspotSourceEntry implements HotspotSourceEntry {

        @NotNull
        private final String gene;
        @Nullable
        private final String transcript;
        @NotNull
        private final String proteinAnnotation;

        public TestHotspotSourceEntry(@NotNull final String gene, @Nullable final String transcript,
                @NotNull final String proteinAnnotation) {
            this.gene = gene;
            this.transcript = transcript;
            this.proteinAnnotation = proteinAnnotation;
        }

        @NotNull
        @Override
        public String gene() {
            return gene;
        }

        @Nullable
        @Override
        public String transcript() {
            return transcript;
        }

        @NotNull
        @Override
        public String proteinAnnotation() {
            return proteinAnnotation;
        }
    }
}