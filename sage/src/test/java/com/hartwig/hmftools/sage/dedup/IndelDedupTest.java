package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_LENGTH;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.addLocalPhaseSet;
import static com.hartwig.hmftools.sage.common.TestUtils.createVariant;
import static com.hartwig.hmftools.sage.common.TestUtils.setTumorQuality;
import static com.hartwig.hmftools.sage.dedup.IndelDeduper.buildAltBasesString;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.vcf.VcfTags.DEDUP_INDEL_FILTER;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.sage.old.IndexedBases;
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
    private static final String CHR_1_REF_BASES = "XAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGTAAAAGGGGCCCCTTTTAAAACCCCGGGGTTTTACGT";

    public IndelDedupTest()
    {
        mRefGenome = new MockRefGenome();
        mIndelDeduper = new IndelDeduper(mRefGenome, DEFAULT_READ_LENGTH);

        mRefGenome.RefGenomeMap.put(CHR_1, CHR_1_REF_BASES);
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

    private static SimpleVariant createSimpleVariant(int position, final String ref, final String alt)
    {
        return new SimpleVariant(CHR_1, position, ref, alt);
    }

    private static SageVariant createSageVariant(
            final String chromosome, int position, int index, final String varReadBases, final String ref, final String alt,
            final int localPhaseSet)
    {
        int leftCoreIndex = index - MIN_CORE_DISTANCE;
        int rightCoreIndex = index + alt.length() - 1 + MIN_CORE_DISTANCE;

        IndexedBases indexBases = new IndexedBases(
                16, index, leftCoreIndex, rightCoreIndex, DEFAULT_READ_CONTEXT_FLANK_SIZE, varReadBases.getBytes());

        SageVariant variant = createVariant(chromosome, position, ref, alt, indexBases);

        addLocalPhaseSet(variant, localPhaseSet, 1);
        setTumorQuality(variant, 5, 1000);

        return variant;
    }

    @Test
    public void testIndels1()
    {
        // only the first INDEL is required to explain the read context, and so the others are de-duped
        String readBases1 = CHR_1_REF_BASES.substring(0, 21) + CHR_1_REF_BASES.substring(30, 50);

        SageVariant var1 = createSageVariant(
                CHR_1, 20, 20, readBases1,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // SNV further up
        String readBases2 = CHR_1_REF_BASES.substring(0, 36) + "A" + CHR_1_REF_BASES.substring(37, 50);

        SageVariant var2 = createSageVariant(
                CHR_1, 36, 36, readBases2,
                CHR_1_REF_BASES.substring(36, 37), "A", 1);

        String readBases3 = CHR_1_REF_BASES.substring(0, 41) + CHR_1_REF_BASES.substring(45, 60);

        SageVariant var3 = createSageVariant(
                CHR_1, 40, 40, readBases3,
                CHR_1_REF_BASES.substring(40, 45), CHR_1_REF_BASES.substring(40, 41), 1);

        mIndelDeduper.dedupVariants(Lists.newArrayList(var1, var2, var3));

        assertTrue(var1.isPassing());
        assertTrue(var2.filters().contains(DEDUP_INDEL_FILTER));
        assertTrue(var3.filters().contains(DEDUP_INDEL_FILTER));
    }

    @Test
    public void testIndels2()
    {
        // all 3 variants are required to explain the read context, and so all remain passing
        String combinedReadBases = CHR_1_REF_BASES.substring(0, 21)
                + CHR_1_REF_BASES.substring(30, 32) + "A" // SNV 1
                + CHR_1_REF_BASES.substring(33, 36) + "A" // SNV 2
                + CHR_1_REF_BASES.substring(37, 41) + CHR_1_REF_BASES.substring(45, 70);

        SageVariant del1 = createSageVariant(
                CHR_1, 20, 20, combinedReadBases,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // germline variant part of the solution but cannot be recovered
        SageVariant var2 = createSageVariant(
                CHR_1, 32, 32, combinedReadBases,
                CHR_1_REF_BASES.substring(32, 33), "A", 1);

        var2.filters().add(MAX_GERMLINE_VAF.filterName());

        // SNV further up
        SageVariant var1 = createSageVariant(
                CHR_1, 36, 36, combinedReadBases,
                CHR_1_REF_BASES.substring(36, 37), "A", 1);

        SageVariant del2 = createSageVariant(
                CHR_1, 40, 40, combinedReadBases,
                CHR_1_REF_BASES.substring(40, 45), CHR_1_REF_BASES.substring(40, 41), 1);

        // initially filtered but still considered
        var1.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());
        del2.filters().add(SoftFilter.FRAGMENT_STRAND_BIAS.filterName());

        // variants not part of the solution
        SageVariant var3 = createSageVariant(
                CHR_1, 38, 38, combinedReadBases,
                CHR_1_REF_BASES.substring(38, 39), "A", 1);

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
        String combinedReadBases = CHR_1_REF_BASES.substring(0, 21) + CHR_1_REF_BASES.substring(30, 41) + CHR_1_REF_BASES.substring(45, 70);

        SageVariant del1 = createSageVariant(
                CHR_1, 20, 20, combinedReadBases,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // SNV further up
        SageVariant var2 = createSageVariant(
                CHR_1, 36, 36, combinedReadBases,
                CHR_1_REF_BASES.substring(36, 37), "A", 1);

        SageVariant del2 = createSageVariant(
                CHR_1, 40, 40, combinedReadBases,
                CHR_1_REF_BASES.substring(40, 45), CHR_1_REF_BASES.substring(40, 41), 1);

        // variants not part of the solution
        SageVariant var4 = createSageVariant(
                CHR_1, 32, 32, combinedReadBases,
                CHR_1_REF_BASES.substring(32, 33), "A", 1);

        SageVariant var5 = createSageVariant(
                CHR_1, 38, 38, combinedReadBases,
                CHR_1_REF_BASES.substring(38, 39), "A", 1);

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

    @Test
    public void testFilteredIndelRecovery()
    {
        // an indel is not recovered when the other variant is largely unphased with it
        String combinedReadBases = CHR_1_REF_BASES.substring(0, 21)
                + CHR_1_REF_BASES.substring(30, 32) + "A" // SNV 1
                + CHR_1_REF_BASES.substring(33, 70);

        String indelBases = CHR_1_REF_BASES.substring(0, 21) + CHR_1_REF_BASES.substring(30, 70);

        SageVariant del = createSageVariant(
                CHR_1, 20, 20, indelBases,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // variant not required by the INDEL
        SageVariant var1 = createSageVariant(
                CHR_1, 32, 32, combinedReadBases,
                CHR_1_REF_BASES.substring(32, 33), "A", 1);

        mIndelDeduper.dedupVariants(Lists.newArrayList(del, var1));

        assertTrue(del.isPassing());
        assertFalse(var1.isPassing());
    }

    @Test
    public void testIndelsLocalPhaseSetConditions()
    {
        // an indel is not recovered when the other variant is largely unphased with it
        String combinedReadBases = CHR_1_REF_BASES.substring(0, 21)
                + CHR_1_REF_BASES.substring(30, 32) + "A" // SNV 1
                + CHR_1_REF_BASES.substring(33, 70);

        SageVariant del = createSageVariant(
                CHR_1, 20, 20, combinedReadBases,
                CHR_1_REF_BASES.substring(20, 30), CHR_1_REF_BASES.substring(20, 21), 1);

        // variant part of the solution
        SageVariant var1 = createSageVariant(
                CHR_1, 32, 32, combinedReadBases,
                CHR_1_REF_BASES.substring(32, 33), "A", 1);

        // initially filtered but recoverable - no longer allow non-passing INDELs to be recovered
        // del.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());

        mIndelDeduper.dedupVariants(Lists.newArrayList(del, var1));

        assertTrue(del.isPassing());
        assertTrue(var1.isPassing());

        // test again but this time with the SNV mostly unphased with the INDEL
        del.filters().clear();
        del.filters().add(SoftFilter.MIN_TUMOR_QUAL.filterName());

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
    }
}
