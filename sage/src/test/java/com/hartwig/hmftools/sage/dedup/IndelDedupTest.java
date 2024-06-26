package com.hartwig.hmftools.sage.dedup;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.addLocalPhaseSet;
import static com.hartwig.hmftools.sage.common.TestUtils.setTumorQuality;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.dedup.IndelDeduper.buildAltBasesString;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantUtils;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;
import com.hartwig.hmftools.sage.filter.SoftFilter;

import org.junit.Test;

public class IndelDedupTest
{
    private final IndelDeduper mIndelDeduper;
    private final MockRefGenome mRefGenome;

    //                                                       10        20        30        40        50        60        70
    //                                             0123456789012345678901234567890123456789012345678901234567890123456789012
    private static final String TEST_REF_BASES = "XAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGTAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGT";

    public IndelDedupTest()
    {
        mRefGenome = new MockRefGenome();
        mIndelDeduper = new IndelDeduper(mRefGenome, DEFAULT_READ_LENGTH);

        mRefGenome.RefGenomeMap.put(CHR_1, TEST_REF_BASES);
    }

    @Test
    public void testVariantAltBuilding()
    {
        List<SimpleVariant> variants = Lists.newArrayList();

        //                 100       110       120       130
        //                 01234567890123456789012345678901
        String refBases = "AAAAGGGGCCCCTTTTAAAACCCCGGGGTTTT";

        variants.add(createSimpleVariant(102, "A", "T"));
        variants.add(createSimpleVariant(126, "GGT", "AAA"));
        variants.add(createSimpleVariant(115, "TAAAAC", "T"));
        variants.add(createSimpleVariant(109, "C", "CGGGG"));
        String altBases = buildAltBasesString(refBases, 100, 131, variants);

        assertEquals(altBases, "AATAGGGGCCGGGGCCTTTTCCCGGAAATTT");
    }

    private static SageVariant createIndelTestVariant(
            int position, int index, final String varReadBases, final String ref, final String alt, final int localPhaseSet)
    {
        SimpleVariant variant = createSimpleVariant(position, ref, alt);

        // for these tests, the core is +/- 2 around the index no matter what the indel or variant
        String leftFlank = varReadBases.substring(index - MIN_CORE_DISTANCE - DEFAULT_FLANK_LENGTH, index - 2);

        int rightCoreStart = index + alt.length();
        // int rightCoreStart = index + max(ref.length(), alt.length());
        int coreEndIndex = rightCoreStart + MIN_CORE_DISTANCE;
        String rightFlank = varReadBases.substring(coreEndIndex, coreEndIndex + DEFAULT_FLANK_LENGTH);

        String leftCore = varReadBases.substring(index - MIN_CORE_DISTANCE, index);
        String rightCore = varReadBases.substring(rightCoreStart, rightCoreStart + MIN_CORE_DISTANCE);

        /*
        String leftFlank = varReadBases.substring(0, DEFAULT_FLANK_LENGTH);
        String rightFlank = varReadBases.substring(varReadBases.length() - DEFAULT_FLANK_LENGTH);
        String leftCore = varReadBases.substring(DEFAULT_FLANK_LENGTH, index);
        int rightCoreStart = index + max(ref.length(), alt.length());
        String rightCore = varReadBases.substring(rightCoreStart, varReadBases.length() - DEFAULT_FLANK_LENGTH);
        */

        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore, leftFlank, rightFlank);

        SageVariant sageVariant = VariantUtils.createSageVariant(readContext);

        addLocalPhaseSet(sageVariant, localPhaseSet, 1);
        setTumorQuality(sageVariant, 5, 1000);

        return sageVariant;
    }

    // CLEAN-UP: need to fix up or change the variant read contexts
    /*
    @Test
    public void testIndels1()
    {
        // only the first INDEL is required to explain the read context, and so the others are de-duped
        String readBases1 = TEST_REF_BASES.substring(0, 21) + TEST_REF_BASES.substring(30, 50);

        SageVariant var1 = createIndelTestVariant(
                20, 20, readBases1,
                TEST_REF_BASES.substring(20, 30), TEST_REF_BASES.substring(20, 21), 1);

        // SNV further up
        String readBases2 = TEST_REF_BASES.substring(0, 36) + "A" + TEST_REF_BASES.substring(37, 50);

        SageVariant var2 = createIndelTestVariant(
                36, 36, readBases2,
                TEST_REF_BASES.substring(36, 37), "A", 1);

        String readBases3 = TEST_REF_BASES.substring(0, 41) + TEST_REF_BASES.substring(45, 60);

        SageVariant var3 = createIndelTestVariant(
                40, 40, readBases3,
                TEST_REF_BASES.substring(40, 45), TEST_REF_BASES.substring(40, 41), 1);

        mIndelDeduper.dedupVariants(Lists.newArrayList(var1, var2, var3));

        assertTrue(var1.isPassing());
        assertTrue(var2.filters().contains(DEDUP_INDEL_FILTER));
        assertTrue(var3.filters().contains(DEDUP_INDEL_FILTER));
    }

    @Test
    public void testIndels2()
    {
        // all 3 variants are required to explain the read context, and so all remain passing
        String combinedReadBases = TEST_REF_BASES.substring(0, 21)
                + TEST_REF_BASES.substring(30, 32) + "A" // SNV 1
                + TEST_REF_BASES.substring(33, 36) + "A" // SNV 2
                + TEST_REF_BASES.substring(37, 41) + TEST_REF_BASES.substring(45, 70);

        SageVariant del1 = createIndelTestVariant(
                20, 20, combinedReadBases,
                TEST_REF_BASES.substring(20, 30), TEST_REF_BASES.substring(20, 21), 1);

        // germline variant part of the solution but cannot be recovered
        SageVariant var2 = createIndelTestVariant(
                32, 32, combinedReadBases,
                TEST_REF_BASES.substring(32, 33), "A", 1);

        var2.filters().add(MAX_GERMLINE_VAF.filterName());

        // SNV further up
        SageVariant var1 = createIndelTestVariant(
                36, 36, combinedReadBases,
                TEST_REF_BASES.substring(36, 37), "A", 1);

        SageVariant del2 = createIndelTestVariant(
                40, 40, combinedReadBases,
                TEST_REF_BASES.substring(40, 45), TEST_REF_BASES.substring(40, 41), 1);

        // initially filtered but still considered
        var1.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());
        del2.filters().add(SoftFilter.FRAGMENT_STRAND_BIAS.filterName());

        // variants not part of the solution
        SageVariant var3 = createIndelTestVariant(
                38, 38, combinedReadBases,
                TEST_REF_BASES.substring(38, 39), "A", 1);

        mIndelDeduper.dedupVariants(Lists.newArrayList(del1, var1, del2, var2, var3));

        assertTrue(del1.isPassing());
        assertTrue(var1.isPassing());
        assertTrue(del2.isPassing());

        assertTrue(var2.filters().contains(MAX_GERMLINE_VAF.filterName()));
        assertFalse(var2.filters().contains(DEDUP_INDEL_FILTER));

        assertTrue(var3.filters().contains(DEDUP_INDEL_FILTER));
    }

    @Test
    public void testIndels3()
    {
        // only the 2 INDELs are required to explain the read context, and so all remain passing
        String combinedReadBases = TEST_REF_BASES.substring(0, 21) + TEST_REF_BASES.substring(30, 41) + TEST_REF_BASES.substring(45, 70);

        SageVariant del1 = createIndelTestVariant(
                20, 20, combinedReadBases,
                TEST_REF_BASES.substring(20, 30), TEST_REF_BASES.substring(20, 21), 1);

        // SNV further up
        SageVariant var2 = createIndelTestVariant(
                36, 36, combinedReadBases,
                TEST_REF_BASES.substring(36, 37), "A", 1);

        SageVariant del2 = createIndelTestVariant(
                40, 40, combinedReadBases,
                TEST_REF_BASES.substring(40, 45), TEST_REF_BASES.substring(40, 41), 1);

        // variants not part of the solution
        SageVariant var4 = createIndelTestVariant(
                32, 32, combinedReadBases,
                TEST_REF_BASES.substring(32, 33), "A", 1);

        SageVariant var5 = createIndelTestVariant(
                38, 38, combinedReadBases,
                TEST_REF_BASES.substring(38, 39), "A", 1);

        // initially filtered but still considered
        del2.filters().add(SoftFilter.FRAGMENT_STRAND_BIAS.filterName());
        var4.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());
        var5.filters().add(SoftFilter.MIN_TUMOR_VAF.filterName());

        mIndelDeduper.dedupVariants(Lists.newArrayList(del1, var2, del2, var4, var5));

        assertTrue(del1.isPassing());
        assertTrue(del2.isPassing());
        assertTrue(var2.filters().contains(DEDUP_INDEL_FILTER));
        assertFalse(var4.isPassing());
        assertFalse(var5.isPassing());
    }
    */

    @Test
    public void testFilteredIndelRecovery()
    {
        // an indel is not recovered when the other variant is largely unphased with it
        String combinedReadBases = TEST_REF_BASES.substring(0, 21)
                + TEST_REF_BASES.substring(30, 32) + "A" // SNV 1
                + TEST_REF_BASES.substring(33, 70);

        String indelBases = TEST_REF_BASES.substring(0, 21) + TEST_REF_BASES.substring(30, 70);

        SageVariant del = createIndelTestVariant(
                20, 20, indelBases,
                TEST_REF_BASES.substring(20, 30), TEST_REF_BASES.substring(20, 21), 1);

        // variant not required by the INDEL
        SageVariant var1 = createIndelTestVariant(
                32, 32, combinedReadBases,
                TEST_REF_BASES.substring(32, 33), "A", 1);

        mIndelDeduper.dedupVariants(Lists.newArrayList(del, var1));

        assertTrue(del.isPassing());
        assertFalse(var1.isPassing());
    }

    @Test
    public void testIndelsLocalPhaseSetConditions()
    {
        // an indel is not recovered when the other variant is largely unphased with it
        String combinedReadBases = TEST_REF_BASES.substring(0, 21)
                + TEST_REF_BASES.substring(30, 32) + "A" // SNV 1
                + TEST_REF_BASES.substring(33, 70);

        SageVariant del = createIndelTestVariant(
                20, 20, combinedReadBases, TEST_REF_BASES.substring(20, 30), TEST_REF_BASES.substring(20, 21), 1);

        // variant part of the solution
        SageVariant var1 = createIndelTestVariant(
                32, 32, combinedReadBases, TEST_REF_BASES.substring(32, 33), "A", 1);

        // initially filtered but recoverable - no longer allow non-passing INDELs to be recovered
        // del.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());

        mIndelDeduper.dedupVariants(Lists.newArrayList(del, var1));

        assertTrue(del.isPassing());
        assertTrue(var1.isPassing());

        // test again but this time with the SNV mostly unphased with the INDEL
        del.filters().clear();
        del.filters().add(SoftFilter.MIN_TUMOR_QUAL);

        var1.tumorReadCounters().get(0).addLocalPhaseSet(2, 100, 0);

        mIndelDeduper.dedupVariants(Lists.newArrayList(del, var1));

        assertFalse(del.isPassing());
        assertTrue(var1.isPassing());
    }

    @Test
    public void testRepeatedInserts()
    {
        // 4 inserts all with the same read context

        //                               10        20
        //                     0123456789012345678901234567890
        String chr2RefBases = "XGATCGATCGATCTCTCTCTCGATCGATCGA";
        mRefGenome.RefGenomeMap.put(CHR_2, chr2RefBases);

        // CLEAN-UP: this sets chromosome 2 to not interfere with the other tests - can this use chr 1 as well??

        /*
        String insertedBases = "CT";

        String readBases = chr2RefBases.substring(9, 13) + insertedBases + chr2RefBases.substring(13, 23);

        int position = 12;
        int posReadIndex = 3;
        SageVariant insert1 = createSageVariant(
                CHR_2, position, posReadIndex, readBases, chr2RefBases.substring(position, position + 1),
                chr2RefBases.substring(position, position + 1) + insertedBases, 1);

        position = 14;
        SageVariant insert2 = createSageVariant(
                CHR_2, position, posReadIndex, readBases, chr2RefBases.substring(position, position + 1),
                chr2RefBases.substring(position, position + 1) + insertedBases, 1);

        position = 16;
        SageVariant insert3 = createSageVariant(
                CHR_2, position, posReadIndex, readBases, chr2RefBases.substring(position, position + 1),
                chr2RefBases.substring(position, position + 1) + insertedBases, 1);

        position = 18;
        SageVariant insert4 = createSageVariant(
                CHR_2, position, posReadIndex, readBases, chr2RefBases.substring(position, position + 1),
                chr2RefBases.substring(position, position + 1) + insertedBases, 1);

        // must be same LPS
        addLocalPhaseSet(insert1, 1, 1);
        addLocalPhaseSet(insert2, 1, 1);
        addLocalPhaseSet(insert3, 1, 1);
        addLocalPhaseSet(insert1, 1, 1);

        List<SageVariant> variants = Lists.newArrayList(insert1, insert2, insert3, insert4);

        // will keep left most since all are of equal length

        mIndelDeduper.dedupVariants(variants);

        assertTrue(insert1.isPassing());
        assertFalse(insert2.isPassing());
        assertFalse(insert3.isPassing());
        assertFalse(insert4.isPassing());
        */
    }
}
