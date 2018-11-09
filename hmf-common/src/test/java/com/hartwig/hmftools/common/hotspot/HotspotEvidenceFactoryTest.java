package com.hartwig.hmftools.common.hotspot;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
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

    private void putBase(@NotNull final String base, int tumorCount, int tumourQuality, int normalCount) {
        this.tumorCount.put(base, tumorCount);
        this.tumorQuality.put(base, tumourQuality);
        this.normalCount.put(base, normalCount);
    }

    private void clearCountsAndQuality() {
        tumorCount.clear();
        tumorQuality.clear();
        normalCount.clear();
    }

    @Test
    public void testRemoveDuplicates() {
        putBase("+GTTT", 7, 17, 5);
        putBase("-GTTT", 7, 17, 5);
        final List<Pileup> tumor = Lists.newArrayList(tumor(10));
        final List<Pileup> normal = Collections.emptyList();

        final ListMultimap<Chromosome, VariantHotspot> hotspots = ArrayListMultimap.create();
        assertEquals(2, new HotspotEvidenceFactory(hotspots).evidence(tumor, normal).size());

        final VariantHotspot insert = create(POS, REF, "GTTT");
        final VariantHotspot delete = create(POS, "GTTT", REF);

        hotspots.put(HumanChromosome.fromString(CHROM), insert);
        assertEquals(2, new HotspotEvidenceFactory(hotspots).evidence(tumor, normal).size());

        hotspots.put(HumanChromosome.fromString(CHROM), insert);
        assertEquals(2, new HotspotEvidenceFactory(hotspots).evidence(tumor, normal).size());

        hotspots.put(HumanChromosome.fromString(CHROM), delete);
        assertEquals(2, new HotspotEvidenceFactory(hotspots).evidence(tumor, normal).size());
    }

    @Test
    public void testInframeIndelCreation() {

        final ListMultimap<Chromosome, VariantHotspot> hotspots = ArrayListMultimap.create();
        final HotspotEvidenceFactory victim = new HotspotEvidenceFactory(hotspots);

        Pileup tumor = tumor(10);
        List<HotspotEvidence> result = victim.evidence(Lists.newArrayList(tumor), Collections.emptyList());
        assertEquals(0, result.size());

        putBase("+GTT", 7, 17, 5);
        result = victim.evidence(Lists.newArrayList(tumor), Collections.emptyList());
        assertEquals(0, result.size());

        putBase("-GTT", 7, 17, 5);
        result = victim.evidence(Lists.newArrayList(tumor), Collections.emptyList());
        assertEquals(0, result.size());

        putBase("+GTTT", 7, 17, 5);
        tumor = tumor(10);
        result = victim.evidence(Lists.newArrayList(tumor), Collections.emptyList());
        assertEquals(1, result.size());
        assertEquals("G", result.get(0).ref());
        assertEquals("GTTT", result.get(0).alt());
        assertEquals(0, result.get(0).normalReads());
        assertEquals(0, result.get(0).normalEvidence());

        putBase("-GACC", 8, 18, 6);
        tumor = tumor(10);
        Pileup normal = normal(5);
        result = victim.evidence(Lists.newArrayList(tumor), Lists.newArrayList(normal));
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

        putBase("A", 4, 14, 2);
        assertHotspotCreating(REF, "A", "A", 10, 32);

        putBase("T", 5, 15, 3);
        assertHotspotCreating(REF, "T", "T", 20, 31);

        putBase("C", 6, 16, 4);
        assertHotspotCreating(REF, "C", "C", 10, 31);

        clearCountsAndQuality();
        putBase("+GTTT", 7, 17, 5);
        assertHotspotCreating(REF, "GTTT", "+GTTT", 10, 31);

        clearCountsAndQuality();
        putBase("-GACC", 8, 18, 6);
        assertHotspotCreating("GACC", REF, "-GACC", 10, 31);
    }

    private void assertHotspotCreating(String ref, String alt, String key, int tumorReads, int normalReads) {
        final ListMultimap<Chromosome, VariantHotspot> hotspots = ArrayListMultimap.create();
        hotspots.put(HumanChromosome.fromString(CHROM), create(POS, ref, alt));

        final HotspotEvidenceFactory factory = new HotspotEvidenceFactory(hotspots);

        final Pileup tumor = tumor(tumorReads);
        final Pileup normal = normal(normalReads);

        final HotspotEvidence victim = factory.evidence(Lists.newArrayList(tumor), Lists.newArrayList(normal)).get(0);
        assertEquals(POS, victim.position());
        assertEquals(ref, victim.ref());
        assertEquals(alt, victim.alt());
        assertEquals((int) tumorQuality.get(key), victim.qualityScore());
        assertEquals((int) tumorCount.get(key), victim.tumorEvidence());
        assertEquals(tumorReads, victim.tumorReads());
        assertEquals((int) normalCount.get(key), victim.normalEvidence());
        assertEquals(normalReads, victim.normalReads());

        final HotspotEvidence victim2 = factory.evidence(Lists.newArrayList(tumor), Collections.emptyList()).get(0);
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
    private VariantHotspot create(long position, @NotNull final String ref, @NotNull final String alt) {
        return ImmutableVariantHotspot.builder().chromosome(CHROM).position(position).ref(ref).alt(alt).build();
    }

    @NotNull
    private Pileup tumor(int readCount) {
        return PileupFileTest.create(CHROM, POS, readCount, REF, tumorCount, tumorQuality);
    }

    @NotNull
    private Pileup normal(int readCount) {
        return PileupFileTest.create(CHROM, POS, readCount, REF, normalCount, Maps.newHashMap());
    }

}
