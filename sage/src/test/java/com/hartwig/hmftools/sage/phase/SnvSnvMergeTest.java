package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.FilterConfig;
import com.hartwig.hmftools.sage.config.ImmutableFilterConfig;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class SnvSnvMergeTest {

    private SnvSnvMerge victim;
    private SageVariantList consumer;
    private byte[] refBases;

    @Before
    public void setup() {
        final FastaSequenceIndex index =
                new FastaSequenceIndex(new File(getClass().getResource("/refsequence/refsequence.fasta.fai").getPath()));
        final IndexedFastaSequenceFile sequence =
                new IndexedFastaSequenceFile(new File(getClass().getResource("/refsequence/refsequence.fasta").getPath()), index);

        final FilterConfig config = ImmutableFilterConfig.builder()
                .hardFilter(false)
                .hardMinTumorAltSupport(0)
                .hardMinTumorQual(0)
                .softHotspotFilter(FilterConfig.NO_FILTER)
                .softPanelFilter(FilterConfig.NO_FILTER)
                .softWideFilter(FilterConfig.NO_FILTER)
                .build();

        consumer = new SageVariantList();
        final MnvFactory mnvFactory =
                new MnvFactory(sequence, new SageVariantFactory(config, Lists.newArrayList(), Lists.newArrayList()));
        victim = new SnvSnvMerge(config, consumer, mnvFactory);
        refBases = sequence.getSubsequenceAt("1", 1, 60).getBases();
    }

    @Test
    public void testSnvsArePhased() {
        byte[] altBases = Arrays.copyOf(refBases, refBases.length);
        altBases[3 - 1] = (byte) 'C';
        altBases[4 - 1] = (byte) 'C';
        victim.accept(createSnv("1", 3, altBases, 1));
        victim.accept(createSnv("1", 4, altBases, 1));
        victim.flush();

        assertEquals(consumer.size(), 3);
        assertVariant(3, "GA", "CC", false, consumer.get(0));
        assertVariant(3, "G", "C", true, consumer.get(1));
        assertVariant(4, "A", "C", true, consumer.get(2));
    }

    @Test
    public void testSnvsAreNotPhased() {
        byte[] altBases = Arrays.copyOf(refBases, refBases.length);
        altBases[3 - 1] = (byte) 'C';
        altBases[4 - 1] = (byte) 'C';
        victim.accept(createSnv("1", 3, altBases, 0));
        victim.accept(createSnv("1", 4, altBases, 0));
        victim.flush();

        assertEquals(consumer.size(), 2);
        assertVariant(3, "G", "C", false, consumer.get(0));
        assertVariant(4, "A", "C", false, consumer.get(1));
    }

    @Test
    public void testSnvsHaveASingleGap() {
        byte[] altBases = Arrays.copyOf(refBases, refBases.length);
        altBases[3 - 1] = (byte) 'C';
        altBases[5 - 1] = (byte) 'C';
        victim.accept(createSnv("1", 3, altBases, 1));
        victim.accept(createSnv("1", 5, altBases, 1));
        victim.flush();

        assertEquals(consumer.size(), 3);
        assertVariant(3, "GAT", "CAC", false, consumer.get(0));
        assertVariant(3, "G", "C", true, consumer.get(1));
        assertVariant(5, "T", "C", true, consumer.get(2));
    }

    @Test
    public void testSnvsHaveASingleGapWithUnphasedVariantInIt() {
        byte[] altBasesPhase1 = Arrays.copyOf(refBases, refBases.length);
        byte[] altBasesPhase2 = Arrays.copyOf(refBases, refBases.length);
        altBasesPhase1[3 - 1] = (byte) 'C';
        altBasesPhase2[4 - 1] = (byte) 'C';
        altBasesPhase1[5 - 1] = (byte) 'C';
        victim.accept(createSnv("1", 3, altBasesPhase1, 1));
        victim.accept(createSnv("1", 4, altBasesPhase2, 0));
        victim.accept(createSnv("1", 5, altBasesPhase1, 1));
        victim.flush();

        assertEquals(consumer.size(), 4);
        assertVariant(3, "GAT", "CAC", false, consumer.get(0));
        assertVariant(3, "G", "C", true, consumer.get(1));
        assertVariant(4, "A", "C", false, consumer.get(2));
        assertVariant(5, "T", "C", true, consumer.get(3));
    }

    @Test
    public void testSnvsHaveALargerGap() {
        byte[] altBases = Arrays.copyOf(refBases, refBases.length);
        altBases[3 - 1] = (byte) 'C';
        altBases[6 - 1] = (byte) 'C';
        victim.accept(createSnv("1", 3, altBases, 1));
        victim.accept(createSnv("1", 6, altBases, 1));
        victim.flush();

        assertEquals(consumer.size(), 2);
        assertVariant(3, "G", "C", false, consumer.get(0));
        assertVariant(6, "A", "C", false, consumer.get(1));
    }

    @Test
    public void testLargeMNV() {
        byte[] altBases = Arrays.copyOf(refBases, refBases.length);
        altBases[3 - 1] = (byte) 'C';
        altBases[4 - 1] = (byte) 'C';
        altBases[5 - 1] = (byte) 'C';
        altBases[7 - 1] = (byte) 'C';
        victim.accept(createSnv("1", 3, altBases, 1));
        victim.accept(createSnv("1", 4, altBases, 1));
        victim.accept(createSnv("1", 5, altBases, 1));
        victim.accept(createSnv("1", 7, altBases, 1));
        victim.flush();

        assertEquals(5, consumer.size());
        assertVariant(3, "GATAA", "CCCAC", false, consumer.get(0));
        assertVariant(3, "G", "C", true, consumer.get(1));
        assertVariant(4, "A", "C",true , consumer.get(2));
        assertVariant(5, "T", "C", true, consumer.get(3));
        assertVariant(7, "A", "C", true, consumer.get(4));
    }

    @Test
    public void testLargeIndelBeforeSnv() {
        byte[] altBases = Arrays.copyOf(refBases, refBases.length);
        altBases[4 - 1] = (byte) 'C';
        altBases[7 - 1] = (byte) 'C';

        victim.accept(DedupSnvTest.createIndel("1", 3, new String(refBases), new String(altBases), 10, 1, 0));
        victim.accept(createSnv("1", 4, altBases, 1));
        victim.accept(createSnv("1", 7, altBases, 1));
        victim.flush();
        assertEquals(3, consumer.size());

    }

    static void assertVariant(long position, String ref, String alt, final boolean filtered, SageVariant variant) {
        assertEquals(position, variant.position());
        assertEquals(ref, variant.normal().ref());
        assertEquals(alt, variant.normal().alt());
        if (filtered) {
            assertFalse(variant.filters().isEmpty());
        } else {
            assertTrue(variant.filters().isEmpty());
        }
    }

    @NotNull
    private SageVariant createSnv(@NotNull final String chromosome, long position, @NotNull final byte[] altBases, int localPhaseSet) {

        int readIndex = (int) position - 1;
        final String ref = new String(refBases, readIndex, 1);
        final String alt = new String(altBases, readIndex, 1);

        VariantHotspot variant = ImmutableVariantHotspotImpl.builder().chromosome(chromosome).position(position).ref(ref).alt(alt).build();
        AltContext altContext = new AltContext("SAMPLE", variant);
        SageVariant result = new SageVariant(SageVariantTier.WIDE, Sets.newHashSet(), altContext, Lists.newArrayList(altContext));
        result.localPhaseSet(localPhaseSet);
        altContext.setPrimaryReadContext(new ReadContextCounter(variant,
                new ReadContext(Strings.EMPTY, 0, readIndex, readIndex, readIndex, 2, altBases)));
        return result;
    }

    static class SageVariantList extends ArrayList<SageVariant> implements Consumer<SageVariant> {
        @Override
        public void accept(final SageVariant sageVariant) {
            this.add(sageVariant);
        }
    }

}
