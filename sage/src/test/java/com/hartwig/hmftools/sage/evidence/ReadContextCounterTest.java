package com.hartwig.hmftools.sage.evidence;

import static org.junit.Assert.assertEquals;

import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;

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

    // convert to using MockRefGenome
    private static final IndexedBases REF_BASES = new IndexedBases(550, 0, "TGTTTCTGTTTC".getBytes());

    private static final QualityCalculator QUALITY_CALCULATOR = new QualityCalculator(CONFIG.Quality, RECALIBRATION, REF_BASES);

    @Test
    public void testInsertInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("GT").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 5, "TGTTTC", Strings.EMPTY);

        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(555, "3S3M", "TGTTTC", "######");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR,1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testDeleteInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("GT").alt("G").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 4, "TGTTC", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(556, "2S3M", "TGTTC", "#####");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR,1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testSnvInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("A").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 2, "CAT", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(555, "2S1M", "CAT", "#####");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testRefInLeftSoftClipDoesNotContributeToDepth()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("A").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 2,"CAT", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        String quals = buildQualString(new int[] {37, 37, 37});

        final SAMRecord record = buildSamRecord(555, "2S1M", "CGT", quals);
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(0, victim.depth());
        assertEquals(0, victim.altSupport());
    }

    public static String buildQualString(final int[] quals)
    {
        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < quals.length; ++i)
        {
            sb.append(phredToFastq(37));
        }

        return sb.toString();
    }

    @Test
    public void testMnvInLeftSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("TCG").alt("ATC").position(552).build();
        final ReadContext readContext = createReadContext(552, 2, 0, 6, "GAAAAAT", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(555, "5S3M", "GAAAAATC", "########");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testInsertInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("G").alt("GT").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 5, "TGTTTC", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(553, "2M4S", "TGTTTC", "######");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testDeleteInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("GT").alt("G").position(554).build();
        final ReadContext readContext = createReadContext(554, 1, 0, 4, "TGTTC", Strings.EMPTY);
        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(553, "2M3S", "TGTTC", "#####");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    @Test
    public void testMnvInRightSoftClip()
    {
        final VariantHotspot hotspot = ImmutableVariantHotspotImpl.builder().chromosome("1").ref("TCG").alt("ATC").position(552).build();

        final ReadContext readContext = createReadContext(552, 2, 0, 6, "GAAAAAT", Strings.EMPTY);

        final ReadContextCounter victim =
                new ReadContextCounter(1, hotspot, readContext, TIER, MAX_COVERAGE, 0);

        final SAMRecord record = buildSamRecord(550, "2M6S", "GAAAAATC", "########");
        victim.processRead(record, CONFIG, QUALITY_CALCULATOR, 1);

        assertEquals(1, victim.depth());
        assertEquals(1, victim.altSupport());
    }

    public static ReadContext createReadContext(
            int refPosition, int readIndex, int leftCentreIndex, int rightCentreIndex, String readBases, String microhomology)
    {
        int adjLeftCentreIndex = Math.max(leftCentreIndex, 0);
        int adjRightCentreIndex = Math.min(rightCentreIndex, readBases.length() - 1);
        boolean incompleteCore = adjLeftCentreIndex != leftCentreIndex || adjRightCentreIndex != rightCentreIndex;

        IndexedBases readBasesIndexed = new IndexedBases(refPosition, readIndex, adjLeftCentreIndex, adjRightCentreIndex, 0, readBases.getBytes());

        return new ReadContext(refPosition, "", 0, microhomology, readBasesIndexed, incompleteCore);
    }

    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
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
