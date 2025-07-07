package com.hartwig.hmftools.sage.dedup;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
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
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.filter.SoftFilter;
import com.hartwig.hmftools.sage.filter.VariantFilters;

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
        mIndelDeduper = new IndelDeduper(mRefGenome, new VariantFilters(TEST_CONFIG));

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

        VariantReadContext readContext = createReadContext(variant, leftCore, rightCore, leftFlank, rightFlank);

        SageVariant sageVariant = VariantUtils.createSageVariant(readContext);

        addLocalPhaseSet(sageVariant, localPhaseSet, 1);
        setTumorQuality(sageVariant, 5, 1000);

        return sageVariant;
    }

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
