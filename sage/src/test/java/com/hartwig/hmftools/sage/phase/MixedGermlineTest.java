package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.config.SoftFilter;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.read.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.read.ReadContextCounterTest;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.VariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class MixedGermlineTest
{
    private static final List<HmfTranscriptRegion> TRANSCRIPTS = HmfGenePanelSupplier.allGeneList37();
    @Test
    public void testNonCodingMnv()
    {
        testMixed(100, false);
    }

    @Test
    public void testCodingMnv()
    {
        testMixed(115251156, true);
        testMixed(115251157, false);
        testMixed(115251158, false);
        testMixed(115251159, true);
    }

    @Test
    public void testPsg1()
    {
        final SageVariant somaticSnv = createGermline("19", 43382367, "G", "A");
        final SageVariant mixedMnv = create("19", 43382367, "GT", "AG");
        final SageVariant germlineSnv = create("19", 43382368, "T", "G");

        process("19", somaticSnv, mixedMnv, germlineSnv);

        assertEquals(1, mixedMnv.mixedGermlineImpact());
        assertEquals(1, somaticSnv.mixedGermlineImpact());
        assertEquals(1, germlineSnv.mixedGermlineImpact());

        assertFalse(germlineSnv.isPassing());
        assertFalse(somaticSnv.isPassing());
        assertTrue(mixedMnv.isPassing());
    }

    private void testMixed(int germlinePosition, boolean mvnPass)
    {
        final SageVariant germlineSnv = createGermline("1", germlinePosition, "A", "G");
        final SageVariant somaticSnv = create("1", germlinePosition + 2, "A", "G");
        final SageVariant mixedMnv = create("1", germlinePosition, "ACA", "GCG");

        process("1", germlineSnv, somaticSnv, mixedMnv);

        assertEquals(1, mixedMnv.mixedGermlineImpact());
        assertEquals(1, somaticSnv.mixedGermlineImpact());
        assertEquals(1, germlineSnv.mixedGermlineImpact());

        assertFalse(germlineSnv.isPassing());
        assertEquals(!mvnPass, somaticSnv.isPassing());
        assertEquals(mvnPass, mixedMnv.isPassing());
    }

    private void process(String chromosome, SageVariant... variants)
    {
        final List<SageVariant> consumer = Lists.newArrayList();

        final Phase victim = new Phase(
                TRANSCRIPTS.stream().filter(x -> x.chromosome().equals(chromosome)).collect(Collectors.toList()),
                consumer::add);

        for(SageVariant variant : variants)
        {
            variant.localPhaseSet(1);
            victim.accept(variant);
        }
        assertEquals(0, consumer.size());
        victim.flush();
        assertEquals(variants.length, consumer.size());
    }

    @NotNull
    private static SageVariant createGermline(String chromosome, int position, @NotNull String ref, @NotNull String alt)
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
        final Candidate candidate = new Candidate(VariantTier.PANEL, variant, counter.ReadContext, 0, 0);
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
