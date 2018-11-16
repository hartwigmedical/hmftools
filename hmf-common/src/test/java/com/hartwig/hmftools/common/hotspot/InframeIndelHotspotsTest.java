package com.hartwig.hmftools.common.hotspot;

import static org.junit.Assert.assertTrue;

import java.util.Set;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class InframeIndelHotspotsTest {

    @Test
    public void testDelLength6() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M6D4M", "GATACA");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertContains(indels, 100, "GNNNNNN", "G");
    }

    @Test
    public void testDelLength3() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M3D4M", "GATACA");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertContains(indels, 100, "GNNN", "G");
    }

    @Test
    public void testDelLength2() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M2D4M", "GATACA");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertTrue(indels.isEmpty());
    }

    @Test
    public void testInsLength6() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M6I4M", "GGATCATTACA");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertContains(indels, 100, "G", "GGATCAT");
    }

    @Test
    public void testInsLength5() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M5I4M", "GGATCATACA");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertTrue(indels.isEmpty());
    }

    @Test
    public void testInsLength3() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M3I4M", "GGATTACA");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertContains(indels, 100, "G", "GGAT");
    }

    @Test
    public void testIndel() {
        final SAMRecord record = VariantHotspotEvidenceFactoryTest.buildSamRecord(100, "1M3I1M3D1M", "GGATAT");
        final Set<VariantHotspot> indels = InframeIndelHotspots.findInframeIndelsWithIncorrectRefs(record);
        assertContains(indels, 100, "G", "GGAT");
        assertContains(indels, 101, "ANNN", "A");
    }

    private void assertContains(@NotNull final Set<VariantHotspot> victim, int pos, @NotNull final String ref, @NotNull final String alt) {
        final VariantHotspot expected = ImmutableVariantHotspot.builder().chromosome("*").position(pos).ref(ref).alt(alt).build();
        assertTrue(victim.contains(expected));
    }

}
