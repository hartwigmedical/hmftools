package com.hartwig.hmftools.sage.read;

import static org.junit.Assert.assertEquals;

import java.util.Collections;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.config.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.variant.VariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextCounterTest
{

    private static final String SAMPLE = "sample";
    private static final int MAX_COVERAGE = 1000;
    private static final VariantTier TIER = VariantTier.PANEL;
    private static final SageConfig CONFIG = new SageConfig();
    private static final QualityRecalibrationMap RECALIBRATION = new QualityRecalibrationMap(Collections.emptyList());

    @Test
    public void testInsertInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("GT").position(554).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 554, 1, 0, 5, 0, "TGTTTC".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "3S3M", "TGTTTC", "######");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testDeleteInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("GT").alt("G").position(554).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 554, 1, 0, 4, 0, "TGTTC".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(556, "2S3M", "TGTTC", "#####");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testSnvInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("A").position(554).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 554, 1, 0, 2, 0, "CAT".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "2S1M", "CAT", "#####");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testRefInLeftSoftClipDoesNotContributeToDepth()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("A").position(554).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 554, 1, 0, 2, 0, "CAT".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "2S1M", "CGT", "#####");
        victim.accept(record, CONFIG, 1);

        assertEquals(0, victim.depth());
        assertEquals(0, victim.altSupport());
    }

    @Test
    public void testMnvInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("TCG").alt("ATC").position(552).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 552, 2, 0, 6, 0, "GAAAAAT".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(555, "5S3M", "GAAAAATC", "########");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testInsertInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("GT").position(554).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 554, 1, 0, 5, 0, "TGTTTC".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(553, "2M4S", "TGTTTC", "######");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testDeleteInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("GT").alt("G").position(554).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 554, 1, 0, 4, 0, "TGTTC".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(553, "2M3S", "TGTTC", "#####");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testMnvInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("TCG").alt("ATC").position(552).build();
        final ReadContext readContext = new ReadContext(Strings.EMPTY, 552, 2, 0, 6, 0, "GAAAAAT".getBytes(), Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(SAMPLE, hotspot, readContext, RECALIBRATION, TIER, MAX_COVERAGE, 0, 50, true);

        final SAMRecord record = buildSamRecord(550, "2M6S", "GAAAAATC", "########");
        victim.accept(record, CONFIG, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @NotNull
    public static ReadContext readContext(int refPosition, int readIndex, int leftCentreIndex, int rightCentreIndex, String bases,
            String microhomology)
    {
        return new ReadContext(Strings.EMPTY,
                refPosition,
                readIndex,
                leftCentreIndex,
                rightCentreIndex,
                0,
                bases.getBytes(),
                microhomology);
    }

    @NotNull
    static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
            @NotNull final String qualities)
    {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);
        return record;
    }
}
