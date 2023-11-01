package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.addLocalPhaseSet;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createVariant;
import static com.hartwig.hmftools.sage.common.TestUtils.setTumorQuality;
import static com.hartwig.hmftools.sage.dedup.DedupIndelOld.dedupIndelsOld;
import static com.hartwig.hmftools.sage.dedup.IndelDeduper.buildAltBasesString;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_INDEL_FILTER;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_INDEL_FILTER_OLD;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.filter.SoftFilter;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class IndelDedupTest
{
    private final IndelDeduper mIndelDeduper;
    private final MockRefGenome mRefGenome;

    //                                                       10        20        30        40        50        60        70
    //                                             0123456789012345678901234567890123456789012345678901234567890123456789012
    private static final String CHR_1_REF_BASES = "XAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGTAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGT";

    public IndelDedupTest()
    {
        mRefGenome = new MockRefGenome();
        mIndelDeduper = new IndelDeduper(mRefGenome);

        mRefGenome.RefGenomeMap.put(CHR_1, CHR_1_REF_BASES);
    }

    @Test
    public void testVariantAltBuilding()
    {
        List<VariantHotspot> variants = Lists.newArrayList();

        //                 100       110       120       130
        //                 01234567890123456789012345678901
        String refBases = "AAAAGGGGCCCCTTTTAAAACCCCGGGGTTTT";

        variants.add(createVariantHotspot(102, "A", "T"));
        variants.add(createVariantHotspot(126, "GGT", "AAA"));
        variants.add(createVariantHotspot(115, "TAAAAC", "T"));
        variants.add(createVariantHotspot(109, "C", "CGGGG"));
        String altBases = buildAltBasesString(refBases, 100, 131, variants);

        assertEquals(altBases, "AATAGGGGCCGGGGCCTTTTCCCGGAAATTT");
    }

    private static VariantHotspot createVariantHotspot(int position, final String ref, final String alt)
    {
        return ImmutableVariantHotspotImpl.builder().chromosome(CHR_1).position(position).ref(ref).alt(alt).build();
    }

    private static SageVariant createSageVariant(
            int position, int index, final String varReadBases, final String ref, final String alt, final int localPhaseSet)
    {
        int leftCoreIndex = index - MIN_CORE_DISTANCE;
        int rightCoreIndex = index + alt.length() - 1 + MIN_CORE_DISTANCE;

        IndexedBases indexBases = new IndexedBases(
                16, index, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, varReadBases.getBytes());

        SageVariant variant = createVariant(position, ref, alt, indexBases);

        addLocalPhaseSet(variant, localPhaseSet, 1);
        setTumorQuality(variant, 5, 1000);

        return variant;
    }

    @Test
    public void testIndelsNew1()
    {
        // only the first INDEL is required to explain the read context, and so the others are de-duped
        String readBases1 = CHR_1_REF_BASES.substring(0, 21) + CHR_1_REF_BASES.substring(30, 50);

        SageVariant var1 = createSageVariant(
                20, 20, readBases1,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // SNV further up
        String readBases2 = CHR_1_REF_BASES.substring(0, 36) + "A" + CHR_1_REF_BASES.substring(37, 50);

        SageVariant var2 = createSageVariant(
                36, 36, readBases2,
                CHR_1_REF_BASES.substring(36, 37), "A", 1);

        String readBases3 = CHR_1_REF_BASES.substring(0, 41) + CHR_1_REF_BASES.substring(45, 60);

        SageVariant var3 = createSageVariant(
                40, 40, readBases3,
                CHR_1_REF_BASES.substring(40, 45), CHR_1_REF_BASES.substring(40, 41), 1);

        mIndelDeduper.dedupVariants(Lists.newArrayList(var1, var2, var3));

        assertTrue(var1.isPassing());
        assertTrue(var2.filters().contains(DEDUP_INDEL_FILTER));
        assertTrue(var3.filters().contains(DEDUP_INDEL_FILTER));
    }

    @Test
    public void testIndelsNew2()
    {
        // all 3 variants are required to explain the read context, and so all remain passing
        String combinedReadBases = CHR_1_REF_BASES.substring(0, 21) + CHR_1_REF_BASES.substring(30, 36) + "A"
                + CHR_1_REF_BASES.substring(37, 41) + CHR_1_REF_BASES.substring(45, 70);

        SageVariant var1 = createSageVariant(
                20, 20, combinedReadBases,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // SNV further up
        SageVariant var2 = createSageVariant(
                36, 36, combinedReadBases,
                CHR_1_REF_BASES.substring(36, 37), "A", 1);

        SageVariant var3 = createSageVariant(
                40, 40, combinedReadBases,
                CHR_1_REF_BASES.substring(40, 45), CHR_1_REF_BASES.substring(40, 41), 1);

        // initially filtered but still considered
        var2.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());
        var3.filters().add(SoftFilter.STRAND_BIAS.filterName());

        mIndelDeduper.dedupVariants(Lists.newArrayList(var1, var2, var3));

        assertTrue(var1.isPassing());
        assertTrue(var2.isPassing());
        assertTrue(var3.isPassing());
    }

    @Test
    public void testIndelsNew3()
    {
        // only the 2 INDELs are required to explain the read context, and so all remain passing
        String combinedReadBases = CHR_1_REF_BASES.substring(0, 21) + CHR_1_REF_BASES.substring(30, 41) + CHR_1_REF_BASES.substring(45, 70);

        SageVariant var1 = createSageVariant(
                20, 20, combinedReadBases,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // SNV further up
        SageVariant var2 = createSageVariant(
                36, 36, combinedReadBases,
                CHR_1_REF_BASES.substring(36, 37), "A", 1);

        SageVariant var3 = createSageVariant(
                40, 40, combinedReadBases,
                CHR_1_REF_BASES.substring(40, 45), CHR_1_REF_BASES.substring(40, 41), 1);

        // initially filtered but still considered
        var3.filters().add(SoftFilter.STRAND_BIAS.filterName());

        mIndelDeduper.dedupVariants(Lists.newArrayList(var1, var2, var3));

        assertTrue(var1.isPassing());
        assertTrue(var2.filters().contains(DEDUP_INDEL_FILTER));
        assertTrue(var3.isPassing());
    }

    @Test
    public void testIndels()
    {
        // var1: GATCGATCGA AACTCTCTC TCTCTCTCTC
        // var2: GATCGATCGA AACTC TCTCTCTCTC
        // ref:  GATCGATCGA AAATC TCTCTCTCTC

        // 0123456789  01  2  34  0123456789
        String flank = generateRandomBases(DEFAULT_READ_CONTEXT_FLANK_SIZE);
        String alt1 = "CTCTC";
        String alt2 = "C";
        String readBases1 = flank + "AA" + alt1 + "TC" + "TCTCTCTCTC";
        String readBases2 = flank + "AA" + alt2 + "TC" + "TCTCTCTCTC";
        int leftCoreIndex = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        int index = leftCoreIndex + MIN_CORE_DISTANCE;
        int rightCoreIndex1 = index + alt1.length() - 1 + MIN_CORE_DISTANCE;
        int rightCoreIndex2 = index + alt2.length() - 1 + MIN_CORE_DISTANCE;

        IndexedBases indexBases1 = new IndexedBases(
                12, index, leftCoreIndex, rightCoreIndex1, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases1.getBytes());

        IndexedBases indexBases2 = new IndexedBases(
                16, index, leftCoreIndex, rightCoreIndex2, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases2.getBytes());

        SageVariant var1 = createVariant(12, "A", alt1, indexBases1);
        SageVariant var2 = createVariant(16, "A", alt2, indexBases2);

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        setTumorQuality(var1, 5, 1000);
        setTumorQuality(var2, 5, 800);

        List<SageVariant> variants = Lists.newArrayList(var1, var2);

        dedupIndelsOld(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        // overlapping variants with different read context
        var1 = createVariant(12, "ATG", "A");
        var2 = createVariant(12, "A", "G");

        setTumorQuality(var1, 5, 1000);
        setTumorQuality(var2, 5, 800);
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        dedupIndelsOld(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        var1 = createVariant(12, "A", "ATG");
        var2 = createVariant(12, "A", "G");

        setTumorQuality(var1, 5, 1000);
        setTumorQuality(var2, 5, 800);
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        dedupIndelsOld(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        // delete containing another variant
        var1 = createVariant(12, "ATG", "A");
        var2 = createVariant(14, "G", "A");

        setTumorQuality(var1, 5, 1000);
        setTumorQuality(var2, 5, 800);
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        dedupIndelsOld(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        var1 = createVariant(12, "ATGAT", "A");
        var2 = createVariant(14, "GA", "AT");

        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        dedupIndelsOld(variants);

        assertFalse(var1.isPassing());
        assertTrue(var2.isPassing()); // MNV has longer read context
    }

    @Test
    public void testMultipleIndelDedup()
    {
        // ref:  GATCGATCGA TCTCTCTCTC GATCGATCGA

        // var1: GATCGATCGA AACTCTCTC TCTCTCTCTC
        // var2: GATCGATCGA AACTC TCTCTCTCTC

        // 0123456789  01  2  34  0123456789
        String leftFlank = generateRandomBases(DEFAULT_READ_CONTEXT_FLANK_SIZE);
        String altCore = "TCTCTCTCTCTCTC";
        String rightFlank = leftFlank;

        String readBases1 = leftFlank + altCore + rightFlank;
        String readBases2 = leftFlank + altCore + rightFlank;
        String readBases3 = leftFlank + altCore + rightFlank;
        String readBases4 = leftFlank + altCore + rightFlank;

        int leftCoreIndex = DEFAULT_READ_CONTEXT_FLANK_SIZE;
        int index = leftCoreIndex + MIN_CORE_DISTANCE;
        int rightCoreIndex = leftCoreIndex + altCore.length() - 1;

        IndexedBases indexBases1 = new IndexedBases(
                12, index, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases1.getBytes());

        IndexedBases indexBases2 = new IndexedBases(
                14, index + 2, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases2.getBytes());

        IndexedBases indexBases3 = new IndexedBases(
                16, index + 4, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases3.getBytes());

        IndexedBases indexBases4 = new IndexedBases(
                18, index + 6, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, readBases4.getBytes());

        SageVariant var1 = createVariant(12, "T", "TCT", indexBases1);
        SageVariant var2 = createVariant(14, "T", "TCT", indexBases2);
        SageVariant var3 = createVariant(16, "T", "TCT", indexBases3);
        SageVariant var4 = createVariant(18, "T", "TCT", indexBases4);

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var3, 1, 1);
        addLocalPhaseSet(var4, 1, 1);

        setTumorQuality(var1, 5, 600);
        setTumorQuality(var2, 5, 700);
        setTumorQuality(var3, 5, 800);
        setTumorQuality(var4, 5, 900);

        List<SageVariant> variants = Lists.newArrayList(var1, var2, var3, var4);

        // will keep left most since all are of equal length
        dedupIndelsOld(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());
        assertFalse(var3.isPassing());
        assertFalse(var4.isPassing());
    }

    @Test
    public void testMultipleIndelDedup2()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 150);

        RegionTaskTester tester = new RegionTaskTester();

        tester.PanelRegions.add(new BaseRegion(1, 1000));

        RegionTask task = tester.createRegionTask(region);

        // bases from 19883401
        String refBases = "ATGCTTATGGCACAAGGGTTTGCCTTCATGCCTGTGTTTACTGCACATAGTGCTTCTTCTTTTTTTTTTTTTTTTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCA"
                + "CATTTCCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTACAGAAATGAGAACTGTGGCAAGCCCTGGAACATCACTGTGGAAGAGCAGTAACATTTATGGAAATGAATTGATAACATTCATTAAGGCTAT";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1100)); // need to cover the ref sequence buffer

        List<SAMRecord> reads = Lists.newArrayList(
                createSamRecord(
                        "READ_01", CHR_1,  21,
                        "GCCTTCATGCCTGTGTTTACTGCACATAGTGCTTCTTCTTTTTTTTTTTTTTTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTACAG",
                        "38M1D57M11I45M"),
                createSamRecord(
                        "READ_02", CHR_1,  25,
                        "TCATGCCTGTGTTTACTGCACATAGTGCTTCTTCTTTTTTTTTTTTTTTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTAAAGAAAT",
                        "34M1D57M11I49M"),
                createSamRecord(
                        "READ_03", CHR_1,  34,
                        "TGTTTACTGCACATAGTGCTTCTTCTTTTTTTTTTTTTTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTACAGAAATGAGAACTGTG",
                        "25M2D56M11I59M"),
                createSamRecord(
                        "READ_04", CHR_1,  59,
                        "CTTTTTTTTTTTTTTTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTACAGAAATGAGAACTGTGGCAAGCCCTGGAACATCACTGTG",
                        "60M11I80M"),
                createSamRecord(
                        "READ_05", CHR_1,  73,
                        "TTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTACAGAAATGAGAACTGTGGCAAGCCCTGGAACATCACTGTGGAAGAGCAGTAACA",
                        "46M11I94M"),
                createSamRecord(
                        "READ_06", CHR_1,  73,
                        "TTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCTTAAGTCTGGCATCATACAGTGAGGTTACAGAAATGAGAACTGTGGCAAGCCCTGGAACATCACTGTGGAAGAGCAGTAACA",
                        "46M11I94M"),
                createSamRecord(
                        "READ_06", CHR_1,  1, // 19883391, with TTTCCTGGGA
                        "TGCTTATGGCACAAGGGTTTGCCTTCTTGCCTGTGTTTACTGCACATAGTGCTTATTATTTTTTTTTTTTTTTTTCATTTACCTCAGTAATTGTAAAAGGGCTTCAAGCCACATTTTCTAAGACACTTCTACAAACATTCT",
                        "57M1D59M25S"));

        reads.get(1).setFirstOfPairFlag(false);
        reads.get(1).setReadNegativeStrandFlag(true);
        reads.get(3).setFirstOfPairFlag(false);
        reads.get(3).setReadNegativeStrandFlag(true);
        reads.get(5).setFirstOfPairFlag(false);
        reads.get(5).setReadNegativeStrandFlag(true);

        tester.TumorSamSlicer.ReadRecords.addAll(reads);
        tester.TumorSamSlicer.ReadRecords.addAll(reads); // to get over qual thresholds

        task.run();

        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 116).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 118 && x.isIndel()).findFirst().orElse(null);
        SageVariant var3 = task.getVariants().stream().filter(x -> x.position() == 118 && !x.isIndel()).findFirst().orElse(null);

        assertNotNull(var1);
        assertNotNull(var2);
        assertNotNull(var3);

        assertTrue(var1.isPassing());
        assertTrue(var2.isPassing());
        assertTrue(var3.filters().contains(DEDUP_INDEL_FILTER_OLD));
    }

}
