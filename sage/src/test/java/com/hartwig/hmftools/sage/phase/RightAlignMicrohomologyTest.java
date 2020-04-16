package com.hartwig.hmftools.sage.phase;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantTier;

import org.junit.Test;

public class RightAlignMicrohomologyTest {

    @Test
    public void testRightAlignDel() {

        final String microhomology = "AG";
        final VariantHotspot leftAligned = BufferedProcessorTest.create(100, "CAGAAACCCATGTATGAAGTACAGTGGA", "C");
        final VariantHotspot rightAligned = RightAlignMicrohomology.rightAlignDel(leftAligned, microhomology);

        assertEquals(102, rightAligned.position());
        assertEquals("GAAACCCATGTATGAAGTACAGTGGAAG", rightAligned.ref());
        assertEquals("G", rightAligned.alt());
    }

    @Test
    public void testKit() {
        final List<HmfTranscriptRegion> transcriptRegions =
                HmfGenePanelSupplier.allGeneList37().stream().filter(x -> x.chromosome().equals("4")).collect(Collectors.toList());

        final  List<SageVariant> collector = Lists.newArrayList();
        final RightAlignMicrohomology victim = new RightAlignMicrohomology(collector::add, transcriptRegions);

        final VariantHotspot leftAligned = ImmutableVariantHotspotImpl.builder()
                .chromosome("4")
                .position(55593579)
                .ref("CAGAAACCCATGTATGAAGTACAGTGGA")
                .alt("C")
                .build();

        ReadContextCounter counter = MixedGermlineTest.dummyCounter(leftAligned, "AG");
        SageVariant variant =
                new SageVariant(SageVariantTier.PANEL, leftAligned, Sets.newHashSet(), Lists.newArrayList(), Lists.newArrayList(counter));

        assertTrue(victim.realign(variant));
    }

}
