package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_CORE;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_FULL;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_PARTIAL;

import static org.junit.Assert.assertEquals;

import static htsjdk.samtools.SAMUtils.phredToFastq;

import java.util.Collections;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.pipeline.RegionTask;
import com.hartwig.hmftools.sage.quality.QualityCalculator;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.common.VariantTier;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

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

    @Test
    public void testBasicVariants()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = "XTTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCAAATT"
                + "TTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCAAATT";
        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        // first a read that establishes the SNV
        String readBases = refBases.substring(1, 51) + "T" + refBases.substring(52, 71);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 1, readBases, "70M");
        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        // a read with a DEL before this variant
        readBases = refBases.substring(1, 31)  + refBases.substring(41, 71);
        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 1, readBases, "30M10D30M");
        tester.TumorSamSlicer.ReadRecords.add(read2);
        tester.TumorSamSlicer.ReadRecords.add(read2);

        // now a read with a DEL so the part of the read that supports the variant is in a soft-clipped section
        readBases = refBases.substring(1, 31) + refBases.substring(41, 51)+ "T" + refBases.substring(52, 56);
        SAMRecord read3 = createSamRecord("READ_03", CHR_1, 1, readBases, "30M15S");
        tester.TumorSamSlicer.ReadRecords.add(read3);

        readBases = refBases.substring(46, 51)+ "T" + refBases.substring(52, 92);
        SAMRecord read4 = createSamRecord("READ_04", CHR_1, 62, readBases, "16S30M");
        tester.TumorSamSlicer.ReadRecords.add(read4);

        task.run();

        TestCase.assertEquals(2, task.getVariants().size());
        SageVariant var = task.getVariants().stream().filter(x -> x.position() == 51).findFirst().orElse(null);
        TestCase.assertNotNull(var);
        TestCase.assertEquals(2, var.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(1, var.tumorReadCounters().get(0).counts()[RC_PARTIAL]);
        TestCase.assertEquals(1, var.tumorReadCounters().get(0).counts()[RC_CORE]);
    }
}
