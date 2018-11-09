package com.hartwig.hmftools.common.hotspot;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.pileup.Pileup;
import com.hartwig.hmftools.common.pileup.PileupFileTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class HotspotEvidenceFactoryTest {

    private static final long POS = 1;
    private static final String REF = "G";
    private static final String CHROM = "1";

    private Map<String, Integer> normalCount;
    private Map<String, Integer> tumorCount;
    private Map<String, Integer> tumorQuality;

    @Before
    public void setup() {
        normalCount = Maps.newHashMap();
        tumorCount = Maps.newHashMap();
        tumorQuality = Maps.newHashMap();

    }

    private void putBase(String base, int tumorCount, int tumourQuality, int normalCount) {
        this.tumorCount.put(base, tumorCount);
        this.tumorQuality.put(base, tumourQuality);
        this.normalCount.put(base, normalCount);
    }

    @Test
    public void testInframeIndelCreation() {
        Pileup tumor = tumor(10);
        List<HotspotEvidence> result = HotspotEvidenceFactory.inframeIndelEvidence(tumor, Optional.empty());
        assertEquals(0, result.size());

        putBase("+GTT", 7, 17, 5);
        result = HotspotEvidenceFactory.inframeIndelEvidence(tumor, Optional.empty());
        assertEquals(0, result.size());

        putBase("-GTT", 7, 17, 5);
        result = HotspotEvidenceFactory.inframeIndelEvidence(tumor, Optional.empty());
        assertEquals(0, result.size());

        putBase("+GTTT", 7, 17, 5);
        tumor = tumor(10);
        result = HotspotEvidenceFactory.inframeIndelEvidence(tumor, Optional.empty());
        assertEquals(1, result.size());
        assertEquals("G", result.get(0).ref());
        assertEquals("GTTT", result.get(0).alt());
        assertEquals(0, result.get(0).normalReads());
        assertEquals(0, result.get(0).normalEvidence());

        putBase("-GACC", 8, 18, 6);
        tumor = tumor(10);
        Pileup normal = normal(5);
        result = HotspotEvidenceFactory.inframeIndelEvidence(tumor, Optional.of(normal));
        assertEquals(2, result.size());
        assertEquals("GACC", result.get(1).ref());
        assertEquals("G", result.get(1).alt());
        assertEquals(10, result.get(1).tumorReads());
        assertEquals(8, result.get(1).tumorEvidence());
        assertEquals(5, result.get(1).normalReads());
        assertEquals(6, result.get(1).normalEvidence());
    }

    @Test
    public void testHotspotCreation() {

        putBase("G", 3, 13, 1);
        putBase("A", 4, 14, 2);
        putBase("T", 5, 15, 3);
        putBase("C", 6, 16, 4);
        putBase("+GTTT", 7, 17, 5);
        putBase("-GACC", 8, 18, 6);

        assertHotspotCreating("G", "A", "A", 10, 31);
        assertHotspotCreating("G", "T", "T", 20, 31);
        assertHotspotCreating("G", "C", "C", 10, 31);
        assertHotspotCreating("G", "GTTT", "+GTTT", 10, 31);
        assertHotspotCreating("GACC", "G", "-GACC", 10, 31);
    }

    private void assertHotspotCreating(String ref, String alt, String key, int tumorReads, int normalReads) {
        final VariantHotspot hotspot = create(POS, ref, alt);
        final Pileup tumor = tumor(tumorReads);
        final Pileup normal = normal(normalReads);

        final HotspotEvidence victim = HotspotEvidenceFactory.fromHotspot(hotspot, tumor, Optional.of(normal));
        assertEquals(POS, victim.position());
        assertEquals(ref, victim.ref());
        assertEquals(alt, victim.alt());
        assertEquals((int) tumorQuality.get(key), victim.qualityScore());
        assertEquals((int) tumorCount.get(key), victim.tumorEvidence());
        assertEquals(tumorReads, victim.tumorReads());
        assertEquals((int) normalCount.get(key), victim.normalEvidence());
        assertEquals(normalReads, victim.normalReads());

        final HotspotEvidence victim2 = HotspotEvidenceFactory.fromHotspot(hotspot, tumor, Optional.empty());
        assertEquals(POS, victim2.position());
        assertEquals(ref, victim2.ref());
        assertEquals(alt, victim2.alt());
        assertEquals((int) tumorQuality.get(key), victim2.qualityScore());
        assertEquals((int) tumorCount.get(key), victim2.tumorEvidence());
        assertEquals(tumorReads, victim2.tumorReads());
        assertEquals(0, victim2.normalEvidence());
        assertEquals(0, victim2.normalReads());
    }

    @NotNull
    VariantHotspot create(long position, @NotNull final String ref, @NotNull final String alt) {
        return ImmutableVariantHotspot.builder().chromosome(CHROM).position(position).ref(ref).alt(alt).build();
    }

    @NotNull
    Pileup tumor(int readCount) {
        return PileupFileTest.create(CHROM, POS, readCount, REF, tumorCount, tumorQuality);
    }

    @NotNull
    Pileup normal(int readCount) {
        return PileupFileTest.create(CHROM, POS, readCount, REF, normalCount, Maps.newHashMap());
    }

}
