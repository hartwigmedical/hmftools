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
    public void testDelWithSoftClipSupport()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = "XTTCTGTAGGTTTCAGATGAAATTTTATCCCCCCCCCCCCCCCCCCCCCCCCCCCTTCACTTCAGCAAATT"
                + "TTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCAAATT";
        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        // a read with a DEL before this variant
        String readBases = refBases.substring(1, 31)  + refBases.substring(41, 51) + "T" + refBases.substring(52, 81);
        // readBases = refBases.substring(1, 31)  + refBases.substring(41, 81);
        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 1, readBases, "30M10D40M");
        tester.TumorSamSlicer.ReadRecords.add(read2);
        tester.TumorSamSlicer.ReadRecords.add(read2);

        // now a read with a DEL so the part of the read that supports the variant is in a soft-clipped section
        readBases = refBases.substring(1, 31) + refBases.substring(41, 51) + "T" + refBases.substring(52, 61);
        SAMRecord read3 = createSamRecord("READ_03", CHR_1, 1, readBases, "30M20S");
        tester.TumorSamSlicer.ReadRecords.add(read3);

        task.run();

        TestCase.assertEquals(2, task.getVariants().size());
        SageVariant var = task.getVariants().stream().filter(x -> x.position() == 51).findFirst().orElse(null);
        TestCase.assertNotNull(var);
        TestCase.assertEquals(2, var.tumorReadCounters().get(0).counts()[RC_FULL]);
    }

    @Test
    public void testBrac2MultiVariants()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 450);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = "XCTAACATACAGTTAGCAGCGACAAAAAAAACTCAGTATCAACAACTACCGGTACAAACCTTTCATTGTAATTTTTCAGTTTTGATAAGTGCTTGTTAGTTTATGGAATCT"
                + "CCATATGTTGAATTTTTGTTTTGTTTTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCAAATTTTTAGATCCAGACTTTCAGCCAT"
                + "CTTGTTCTGAGGTGGACCTAATAGGATTTGTCGTTTCTGTTGTGAAAAAAACAGGTAATGCACAATATAGTTAATTTTTTTTATTGATTCTTTTAAAAAACATTGTCTTTTAAAATCT"
                + "CTTATGATTAGTTGGAGCTACCAGTTGGCAAATTTGCTAGCTAACTAGTGATCTGAAAGTAAGCCTCTTTGAACCTCTGATTTTTCATGAAAAGCAATTCTCTC";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        SAMRecord read0 = createSamRecord("READ_00", CHR_1, 171, // 36965
                "TACCAGCCACGGGAGCCCCTTCACTTCAGCAAATTTTTAGATCCAAATAGGATCTAAACAAATAGGATTTGTTTCTGTTGTGAAAAAAACAGGTAATGCACAATATAGT"
                        + "TAATTTTTTTTATTGATTCTTTTAAAAAACATTGTCTTTTAA",
                "47M23D1M1I6M5I10M3D81M");

        tester.TumorSamSlicer.ReadRecords.add(read0);
        tester.TumorSamSlicer.ReadRecords.add(read0);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 87, // 8187
                "AGTGCTTGTTAGTTTATGGAATCTCCATATGTTGAATTTTTGTTTTGTTTTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTT"
                        + "CAGCAAATTTTTAGATCCAAATAGGATCTAAACAAATAGGA",
                "129M22S");

        tester.TumorSamSlicer.ReadRecords.add(read1);

        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 92,
                "TTGTTAGTTTATGGAATCTCCATATGTTGAATTTTTGTTTTGTTTTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCA"
                        + "AATTTTTAGATCCAAATAGGATCTAAACAAATAGGATTTGT",
                "124M27S");

        SAMRecord read3 = createSamRecord("READ_03", CHR_1, 95,
                "TTAGTTTATGGAATCTCCATATGTTGAATTTTTGTTTTGTTTTCTGTAGGTTTCAGATGAAATTTTATTTCAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCAAAT"
                        + "TTTTAGATCCAAATAGGATCTAAACAAATAGGATTTGTTTC",
                "121M30S");

        SAMRecord read4 = createSamRecord("READ_04", CHR_1, 141,
                "TAGGTTTCAGATGAAATTTGATTACAGATTTACCAGCCACGGGAGCCCCTTCACTTCAGCAAATTTTTAGATCCAAATAGGATCTAAACAAATAGGATTTGTTTCTGTTG"
                        + "TGAAAAAAACAGGTAATGCACAATATAGTTAATTTTTTTTA",
                "77M23D1M1I6M5I10M3D51M");

        tester.TumorSamSlicer.ReadRecords.add(read2);
        tester.TumorSamSlicer.ReadRecords.add(read3);
        tester.TumorSamSlicer.ReadRecords.add(read4);

        task.run();

        TestCase.assertEquals(6, task.getVariants().size());
        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 216).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 217).findFirst().orElse(null);
        SageVariant var3 = task.getVariants().stream().filter(x -> x.position() == 241).findFirst().orElse(null);
        SageVariant var4 = task.getVariants().stream().filter(x -> x.position() == 245).findFirst().orElse(null);
        SageVariant var5 = task.getVariants().stream().filter(x -> x.position() == 247).findFirst().orElse(null);
        SageVariant var6 = task.getVariants().stream().filter(x -> x.position() == 257).findFirst().orElse(null);
        TestCase.assertNotNull(var1);
        TestCase.assertNotNull(var2);
        TestCase.assertNotNull(var3);
        TestCase.assertNotNull(var4);
        TestCase.assertNotNull(var5);
        TestCase.assertEquals(5, var1.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(1, var1.tumorReadCounters().get(0).counts()[RC_PARTIAL]);

        TestCase.assertEquals(5, var2.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(1, var2.tumorReadCounters().get(0).counts()[RC_PARTIAL]);

        TestCase.assertEquals(5, var3.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(1, var3.tumorReadCounters().get(0).counts()[RC_PARTIAL]);

        TestCase.assertEquals(4, var4.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(1, var4.tumorReadCounters().get(0).counts()[RC_PARTIAL]);

        TestCase.assertEquals(5, var5.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(1, var5.tumorReadCounters().get(0).counts()[RC_PARTIAL]);

        TestCase.assertEquals(3, var6.tumorReadCounters().get(0).counts()[RC_FULL]);
        TestCase.assertEquals(0, var6.tumorReadCounters().get(0).counts()[RC_PARTIAL]);
    }

}
