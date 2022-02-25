package com.hartwig.hmftools.sage.phase;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.DEFAULT_READ_CONTEXT_FLANK_SIZE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.common.TestUtils.addLocalPhaseSet;
import static com.hartwig.hmftools.sage.common.TestUtils.clearFilters;
import static com.hartwig.hmftools.sage.common.TestUtils.createVariant;
import static com.hartwig.hmftools.sage.common.TestUtils.createVariantHotspot;
import static com.hartwig.hmftools.sage.common.TestUtils.setTumorQuality;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_GERMLINE_DEPTH;
import static com.hartwig.hmftools.sage.dedup.DedupIndel.dedupIndels;
import static com.hartwig.hmftools.sage.dedup.VariantDeduper.longerContainsShorter;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.dedup.DedupMixedGermlineSomatic;
import com.hartwig.hmftools.sage.dedup.DedupSnvMnv;

import org.junit.Test;

public class VariantDedupTest
{
    @Test
    public void testLongerContainsShorter()
    {
        final VariantHotspot mnv = createVariantHotspot(100, "CAC", "TGT");

        assertTrue(longerContainsShorter(createVariantHotspot(100, "C", "T"), mnv));
        assertFalse(longerContainsShorter(createVariantHotspot(100, "C", "C"), mnv));

        assertTrue(longerContainsShorter(createVariantHotspot(101, "A", "G"), mnv));
        assertTrue(longerContainsShorter(createVariantHotspot(102, "C", "T"), mnv));

        assertTrue(longerContainsShorter(createVariantHotspot(100, "CA", "TG"), mnv));
        assertTrue(longerContainsShorter(createVariantHotspot(101, "AC", "GT"), mnv));

        assertTrue(longerContainsShorter(createVariantHotspot(100, "CAC", "TGT"), mnv));
    }

    @Test
    public void testMnvDedup()
    {
        // must overlap, be phased and passing

        // ATG -> CGA
        SageVariant var1 = createVariant(10, "AT", "CG");
        SageVariant var2 = createVariant(11, "TG", "GA");

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
        var1 = createVariant(10, "ATG", "CTA");
        var2 = createVariant(12, "GC", "AT");
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        DedupSnvMnv.dedupMnvOverlaps(variants);

        assertFalse(var1.isPassing());
        assertTrue(var2.isPassing());

        clearFilters(variants);

        // dedup longer assuming all bases have changed
        var1 = createVariant(10, "AAG", "CTA");
        var2 = createVariant(12, "GC", "AT");
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
        SageVariant var1 = createVariant(10, "ATG", "CGA");
        SageVariant var2 = createVariant(10, "A", "C");
        SageVariant var3 = createVariant(11, "T", "G");
        SageVariant var4 = createVariant(12, "G", "A");

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
        var1 = createVariant(10, "A", "G"); // no match
        var2 = createVariant(10, "ATG", "CGA");
        SageVariant var7 = createVariant(11, "TGA", "CCC");
        var3 = createVariant(11, "T", "A");
        var4 = createVariant(11, "T", "G"); // dedup with first MNV
        SageVariant var5 = createVariant(12, "G", "T");
        SageVariant var6 = createVariant(12, "G", "A"); // deup with first MNV
        SageVariant var8 = createVariant(13, "A", "C");

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
        SageVariant var1 = createVariant(10, "ATG", "CGA");
        SageVariant var2 = createVariant(10, "A", "C");
        SageVariant var3 = createVariant(12, "G", "A"); // marked as germline

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var3, 1, 1);

        var3.filters().add(MIN_GERMLINE_DEPTH.filterName());

        List<SageVariant> variants = Lists.newArrayList(var1, var2, var3);

        List<TranscriptData> transcripts = Lists.newArrayList();

        transcripts.add(GeneTestUtils.createTransExons(
                GENE_ID_1, TRANS_ID_1, POS_STRAND, new int[] {7, 20}, 9, 7, 25, true, ""));

        DedupMixedGermlineSomatic deduper = new DedupMixedGermlineSomatic(transcripts);
        deduper.dedupVariants(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());
        assertFalse(var3.isPassing());

        // test again with the MNV spanning a codon so the SNV is kept
        var1 = createVariant(11, "ATG", "CGA");
        var2 = createVariant(12, "T", "G"); // marked as germline
        var3 = createVariant(13, "G", "A");

        // must be same LPS
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);
        addLocalPhaseSet(var3, 1, 1);

        var2.filters().add(MIN_GERMLINE_DEPTH.filterName());

        variants = Lists.newArrayList(var1, var2, var3);

        deduper.dedupVariants(variants);

        assertFalse(var1.isPassing());
        assertFalse(var2.isPassing());
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

        dedupIndels(variants);

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

        dedupIndels(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        var1 = createVariant(12, "A", "ATG");
        var2 = createVariant(12, "A", "G");

        setTumorQuality(var1, 5, 1000);
        setTumorQuality(var2, 5, 800);
        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        dedupIndels(variants);

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

        dedupIndels(variants);

        assertTrue(var1.isPassing());
        assertFalse(var2.isPassing());

        var1 = createVariant(12, "ATGAT", "A");
        var2 = createVariant(14, "GA", "AT");

        addLocalPhaseSet(var1, 1, 1);
        addLocalPhaseSet(var2, 1, 1);

        variants = Lists.newArrayList(var1, var2);

        dedupIndels(variants);

        assertFalse(var1.isPassing());
        assertTrue(var2.isPassing()); // MNV has longer read context
    }

}