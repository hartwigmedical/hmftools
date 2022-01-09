package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.evidence.ReadContextCounter;
import com.hartwig.hmftools.sage.evidence.ReadContextCounterTest;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.VariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MixedGermlineTest
{
    private final List<TranscriptData> mTranscripts;

    public MixedGermlineTest()
    {
        mTranscripts = Lists.newArrayList();

        int[] exonStarts = new int[] {100, 200, 300};
        mTranscripts.add(createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, exonStarts, 110, 100, 310, true, ""));
    }

    @Test
    public void testNonCodingMnv()
    {
        testMixed(50, false);
    }

    @Test
    public void testCodingMnv()
    {
        testMixed(103, true); // start of a codon will pass
        testMixed(104, false);
        testMixed(105, false);
        testMixed(106, true);
    }

    @Test
    public void testPsg1()
    {
        // phases 2 and 0
        final SageVariant somaticSnv = createGermline(CHR_1, 104, "G", "A"); // was 43382367
        final SageVariant mixedMnv = create(CHR_1, 104, "GT", "AG"); // was 43382367
        final SageVariant germlineSnv = create(CHR_1, 105, "T", "G"); // was 43382368

        process(somaticSnv, mixedMnv, germlineSnv);

        assertEquals(1, mixedMnv.mixedGermlineImpact());
        assertEquals(1, somaticSnv.mixedGermlineImpact());
        assertEquals(1, germlineSnv.mixedGermlineImpact());

        assertFalse(germlineSnv.isPassing());
        assertFalse(somaticSnv.isPassing());
        assertTrue(mixedMnv.isPassing());
    }

    private void testMixed(int germlinePosition, boolean mvnPass)
    {
        final SageVariant germlineSnv = createGermline(CHR_1, germlinePosition, "A", "G");
        final SageVariant somaticSnv = create(CHR_1, germlinePosition + 2, "A", "G");
        final SageVariant mixedMnv = create(CHR_1, germlinePosition, "ACA", "GCG");

        process(germlineSnv, somaticSnv, mixedMnv);

        assertEquals(1, mixedMnv.mixedGermlineImpact());
        assertEquals(1, somaticSnv.mixedGermlineImpact());
        assertEquals(1, germlineSnv.mixedGermlineImpact());

        assertFalse(germlineSnv.isPassing());
        assertEquals(!mvnPass, somaticSnv.isPassing());
        assertEquals(mvnPass, mixedMnv.isPassing());
    }

    private void process(SageVariant... variants)
    {
        final List<SageVariant> consumer = Lists.newArrayList();

        final Phase victim = new Phase(mTranscripts, consumer::add);

        for(SageVariant variant : variants)
        {
            variant.localPhaseSet(1);
            victim.accept(variant);
        }
        assertEquals(0, consumer.size());
        victim.flush();
        assertEquals(variants.length, consumer.size());
    }

    private static SageVariant createGermline(final String chromosome, int position, final String ref, final String alt)
    {
        SageVariant result = create(chromosome, position, ref, alt);
        result.filters().add(SoftFilter.MAX_GERMLINE_ALT_SUPPORT.toString());
        return result;
    }

    @NotNull
    private static SageVariant create(String chromosome, int position, @NotNull String ref, @NotNull String alt)
    {
        VariantHotspot variant = ImmutableVariantHotspotImpl.builder().chromosome(chromosome).ref(ref).alt(alt).position(position).build();
        ReadContextCounter counter = dummyCounter(variant, Strings.EMPTY);
        final Candidate candidate = new Candidate(VariantTier.PANEL, variant, counter.readContext(), 0, 0);
        return new SageVariant(candidate, Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(counter));
    }

    @NotNull
    static ReadContextCounter dummyCounter(@NotNull VariantHotspot variant, @NotNull final String microhomology)
    {
        ReadContext dummyReadContext = ReadContextCounterTest.readContext(100, 0, 0, 0, "AAA", microhomology);
        return new ReadContextCounter("SAMPLE",
                variant,
                dummyReadContext,
                new QualityRecalibrationMap(Collections.emptyList()),
                VariantTier.PANEL,
                1000,
                0,
                50,
                false);
    }
}
