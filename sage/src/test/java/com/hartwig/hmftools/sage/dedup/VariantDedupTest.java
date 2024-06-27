package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.sage.common.TestUtils.TEST_CONFIG;
import static com.hartwig.hmftools.sage.common.TestUtils.addLocalPhaseSet;
import static com.hartwig.hmftools.sage.common.TestUtils.clearFilters;
import static com.hartwig.hmftools.sage.common.TestUtils.setTumorQuality;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSageVariant;
import static com.hartwig.hmftools.sage.common.VariantUtils.createSimpleVariant;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_GERMLINE_DEPTH;
import static com.hartwig.hmftools.sage.dedup.VariantDeduper.longerContainsShorter;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SimpleVariant;

import org.junit.Test;

public class VariantDedupTest
{
    @Test
    public void testLongerContainsShorter()
    {
        final SimpleVariant mnv = createSimpleVariant(100, "CAC", "TGT");

        assertTrue(longerContainsShorter(createSimpleVariant(100, "C", "T"), mnv));
        assertFalse(longerContainsShorter(createSimpleVariant(100, "C", "C"), mnv));

        assertTrue(longerContainsShorter(createSimpleVariant(101, "A", "G"), mnv));
        assertTrue(longerContainsShorter(createSimpleVariant(102, "C", "T"), mnv));

        assertTrue(longerContainsShorter(createSimpleVariant(100, "CA", "TG"), mnv));
        assertTrue(longerContainsShorter(createSimpleVariant(101, "AC", "GT"), mnv));

        assertTrue(longerContainsShorter(createSimpleVariant(100, "CAC", "TGT"), mnv));
    }

    @Test
    public void testMnvDedup()
    {
        // must overlap, be phased and passing

        // ATG -> CGA
        SageVariant var1 = createSageVariant(10, "AT", "CG");
        SageVariant var2 = createSageVariant(11, "TG", "GA");

        setTumorQuality(var1, 5, 1000);
        setTumorQuality(var2, 4, 900);

        List<SageVariant> variants = Lists.newArrayList(var1, var2);

        DedupSnvMnv.dedupMnvOverlaps(variants);

        assertTrue(var1.isPassing());
        assertTrue(var2.isPassing());

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        DedupSnvMnv.dedupMnvOverlaps(variants);

        // determined by tumor qual
        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        clearFilters(variants);

        setTumorQuality(var1, 5, 500);

        DedupSnvMnv.dedupMnvOverlaps(variants);

        // determined by tumor qual
        assertFalse(var1.isPassing());
        assertTrue(var2.isPassing());

        clearFilters(variants);

        // dedup shorter if has same number of changed bases
        var1 = createSageVariant(10, "ATG", "CTA");
        var2 = createSageVariant(12, "GC", "AT");
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        DedupSnvMnv.dedupMnvOverlaps(variants);

        assertFalse(var1.isPassing());
        assertTrue(var2.isPassing());

        clearFilters(variants);

        // dedup longer assuming all bases have changed
        var1 = createSageVariant(10, "AAG", "CTA");
        var2 = createSageVariant(12, "GC", "AT");
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        DedupSnvMnv.dedupMnvOverlaps(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());
    }

    @Test
    public void testMnvSnvDedup()
    {
        // ATG -> CGA
        SageVariant var1 = createSageVariant(10, "ATG", "CGA");
        SageVariant var2 = createSageVariant(10, "A", "C");
        SageVariant var3 = createSageVariant(11, "T", "G");
        SageVariant var4 = createSageVariant(12, "G", "A");

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var3, 1, 1);
        addLocalPhaseSet(var4, 1, 1);

        List<SageVariant> variants = Lists.newArrayList(var1, var2, var3, var4);

        DedupSnvMnv.dedupMnvSnvs(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());
        assertFalse(var3.isPassing());
        assertFalse(var4.isPassing());

        // more complicated example with variants interspersed
        var1 = createSageVariant(10, "A", "G"); // no match
        var2 = createSageVariant(10, "ATG", "CGA");
        SageVariant var7 = createSageVariant(11, "TGA", "CCC");
        var3 = createSageVariant(11, "T", "A");
        var4 = createSageVariant(11, "T", "G"); // dedup with first MNV
        SageVariant var5 = createSageVariant(12, "G", "T");
        SageVariant var6 = createSageVariant(12, "G", "A"); // deup with first MNV
        SageVariant var8 = createSageVariant(13, "A", "C");

        // must be same LPS
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var4, 1, 1);
        addLocalPhaseSet(var6, 1, 1);
        addLocalPhaseSet(var7, 2, 1);
        addLocalPhaseSet(var8, 2, 1);

        variants = Lists.newArrayList(var1, var2, var3, var4, var5, var6, var7, var8);

        DedupSnvMnv.dedupMnvSnvs(variants);

        assertTrue(var1.isPassing());
        assertTrue(var2.isPassing());
        assertTrue(var3.isPassing());
        assertFalse(var4.isPassing());
        assertTrue(var5.isPassing());
        assertFalse(var6.isPassing());
        assertTrue(var7.isPassing());
        assertFalse(var8.isPassing());
    }

    @Test
    public void testMixedGermline()
    {
        // ATG -> CGA
        SageVariant var1 = createSageVariant(10, "ATG", "CGA");
        SageVariant var2 = createSageVariant(10, "A", "C");
        SageVariant var3 = createSageVariant(12, "G", "A"); // marked as germline

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var3, 1, 1);

        var3.filters().add(MIN_GERMLINE_DEPTH);

        List<SageVariant> variants = Lists.newArrayList(var1, var2, var3);

        List<TranscriptData> transcripts = Lists.newArrayList();

        transcripts.add(GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {7, 20}, 9, 7, 25, true, ""));

        DedupMixedGermlineSomatic deduper = new DedupMixedGermlineSomatic(transcripts, TEST_CONFIG.Filter);
        deduper.dedupVariants(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());
        assertFalse(var3.isPassing());

        // test again with the MNV spanning a codon so the SNV is kept
        var1 = createSageVariant(11, "ATG", "CGA");
        var2 = createSageVariant(12, "T", "G"); // marked as germline
        var3 = createSageVariant(13, "G", "A");

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var3, 1, 1);

        var2.filters().add(MIN_GERMLINE_DEPTH);

        variants = Lists.newArrayList(var1, var2, var3);

        deduper.dedupVariants(variants);

        assertFalse(var1.isPassing());
        assertFalse(var2.isPassing());
        assertTrue(var3.isPassing());
    }
}