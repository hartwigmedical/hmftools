package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class SnvIndelMergeTest {

    private static final String REF_BASES = "AAGATAAACCATGCAATTTGATATACGGACTTACACAACCATGATGCATTGATCGGACCC";

    private SnvIndelMerge victim;
    private SnvSnvMergeTest.SageVariantList consumer;

    @Before
    public void setup() {
        consumer = new SnvSnvMergeTest.SageVariantList();
        victim = new SnvIndelMerge(consumer);
    }

    @Test
    public void testSnvThenInsert() {
        String altBases = "AAGATAACGTCCCATGCAATTTGATATACGGACTTACACAACCATGATGCATTGATCGGACCC";

        victim.accept(createSnv("1", 8, altBases, 1));
        victim.accept(createIndel("1", 8, altBases, 1, 4, 1));
        victim.flush();

        assertEquals(2, consumer.size());
        SnvSnvMergeTest.assertVariant(8, "A", "C", true, consumer.get(0));
        SnvSnvMergeTest.assertVariant(8, "A", "CGTC", false, consumer.get(1));
    }

    @Test
    public void testInsertThenSnv() {
        String altBases = "AAGATAACGTCCCATGCAATTTGATATACGGACTTACACAACCATGATGCATTGATCGGACCC";

        victim.accept(createIndel("1", 8, altBases, 1, 4, 1));
        victim.accept(createSnv("1", 8, altBases, 1));
        victim.flush();

        assertEquals(2, consumer.size());
        SnvSnvMergeTest.assertVariant(8, "A", "CGTC", false, consumer.get(0));
        SnvSnvMergeTest.assertVariant(8, "A", "C", true, consumer.get(1));
    }

    @Test
    public void testUnphasedSnvThenInsert() {
        String altBases = "AAGATAACGTCCCATGCAATTTGATATACGGACTTACACAACCATGATGCATTGATCGGACCC";

        victim.accept(createSnv("1", 8, altBases, 0));
        victim.accept(createIndel("1", 8, altBases, 1, 4, 0));
        victim.flush();

        assertEquals(2, consumer.size());
        SnvSnvMergeTest.assertVariant(8, "A", "C", false, consumer.get(0));
        SnvSnvMergeTest.assertVariant(8, "A", "CGTC", false, consumer.get(1));
    }

    @Test
    public void testSnvThenDelete() {
        String altBases = "AAGATAACATGCAATTTGATATACGGACTTACACAACCATGATGCATTGATCGGACCC";

        victim.accept(createSnv("1", 8, altBases, 1));
        victim.accept(createIndel("1", 8, altBases, 3, 1, 1));
        victim.flush();

        assertEquals(2, consumer.size());
        SnvSnvMergeTest.assertVariant(8, "A", "C", true, consumer.get(0));
        SnvSnvMergeTest.assertVariant(8, "ACC", "C", false, consumer.get(1));
    }

    @Test
    public void testDeleteThenSnv() {
        String altBases = "AAGATAACATGCAATTTGATATACGGACTTACACAACCATGATGCATTGATCGGACCC";

        victim.accept(createIndel("1", 8, altBases, 3, 1, 1));
        victim.accept(createSnv("1", 8, altBases, 1));
        victim.flush();

        assertEquals(2, consumer.size());
        SnvSnvMergeTest.assertVariant(8, "ACC", "C", false, consumer.get(0));
        SnvSnvMergeTest.assertVariant(8, "A", "C", true, consumer.get(1));
    }

    @NotNull
    private SageVariant createSnv(@NotNull final String chromosome, long position, @NotNull final String altBases, int localPhaseSet) {
        return createIndel(chromosome, position, altBases, 1, 1, localPhaseSet);
    }

    @NotNull
    private SageVariant createIndel(@NotNull final String chromosome, long position, @NotNull final String altBases, int refLength,
            int altLength, int localPhaseSet) {

        int readIndex = (int) position - 1;
        final String ref = REF_BASES.substring(readIndex, readIndex + refLength);
        final String alt = altBases.substring(readIndex, readIndex + altLength);

        VariantHotspot variant = ImmutableVariantHotspotImpl.builder().chromosome(chromosome).position(position).ref(ref).alt(alt).build();
        AltContext altContext = new AltContext("SAMPLE", variant);
        SageVariant result = new SageVariant(SageVariantTier.WIDE, Sets.newHashSet(), altContext, Lists.newArrayList(altContext));
        result.localPhaseSet(localPhaseSet);
        altContext.setPrimaryReadContext(new ReadContextCounter(variant,
                new ReadContext(Strings.EMPTY, 0, readIndex, readIndex, readIndex, 2, altBases.getBytes())));
        return result;
    }

}
