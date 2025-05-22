package com.hartwig.hmftools.lilac.fragment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.NANOS_IN_MILLISECOND;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.mergeFragments;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildSamRecord;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createReadRecord;
import static com.hartwig.hmftools.lilac.read.Read.createRead;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.lilac.read.Read;

import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class FragmentsTest
{
    private static final String TEST_READ_BASES = MockRefGenome.generateRandomBases(100);

    @Test
    public void testIndices()
    {
        assertRange(0, -1, calcAminoAcidIndices(0, 0));
        assertRange(0, -1, calcAminoAcidIndices(0, 1));
        assertRange(0, 0, calcAminoAcidIndices(0, 2));
        assertRange(1, 0, calcAminoAcidIndices(1, 2));
        assertRange(1, 0, calcAminoAcidIndices(2, 2));

        assertRange(1, 1, calcAminoAcidIndices(3, 5));
        assertRange(1, 1, calcAminoAcidIndices(3, 6));
        assertRange(1, 1, calcAminoAcidIndices(3, 7));
        assertRange(1, 2, calcAminoAcidIndices(3, 8));
        assertRange(1, 2, calcAminoAcidIndices(3, 9));
        assertRange(1, 2, calcAminoAcidIndices(3, 10));
    }

    @Test
    public void testAdapterSoftClippedFragments()
    {
        BaseRegion codingRegion = new BaseRegion(100, 1000);

        // firstly a record with soft-clips at both ends which will be used
        SAMRecord record = buildSamRecord(200, "10S60M10S", TEST_READ_BASES.substring(0, 80), "");
        record.setInferredInsertSize(200);

        Read read = createRead(codingRegion, record, true, true);
        assertNotNull(read);

        assertEquals(190, read.PositionStart);
        assertEquals(269, read.PositionEnd);
        assertEquals(10, read.SoftClippedStart);
        assertEquals(10, read.SoftClippedEnd);

        // now restricted at the lower 3' end
        record.setInferredInsertSize(70);
        record.setReadNegativeStrandFlag(true);

        read = createRead(codingRegion, record, true, true);
        assertEquals(200, read.PositionStart);
        assertEquals(269, read.PositionEnd);
        assertEquals(0, read.SoftClippedStart);
        assertEquals(10, read.SoftClippedEnd);

        record.setReadNegativeStrandFlag(false);

        read = createRead(codingRegion, record, true, true);
        assertEquals(190, read.PositionStart);
        assertEquals(259, read.PositionEnd);
        assertEquals(10, read.SoftClippedStart);
        assertEquals(0, read.SoftClippedEnd);
    }

    @Test
    public void testFragmentAddNucleotide()
    {
        List<Integer> indices = Lists.newArrayList(1, 2, 3, 6, 7, 8);
        List<Integer> qualities = Lists.newArrayList(37, 25, 37, 37, 37, 25);
        List<String> nucleotides = Lists.newArrayList("A", "G", "T", "C", "A", "G");

        Fragment fragment = new Fragment(createReadRecord("01"), HLA_A, Sets.newHashSet(HLA_A), indices, qualities, nucleotides);

        fragment.removeLowQualBases();

        assertFalse(fragment.containsNucleotide(2));

        fragment.addNucleotideInfo(5, "G", 30);

        assertTrue(fragment.containsNucleotide(5));
        assertTrue(FragmentUtils.validateLociBases(fragment.id(), fragment.rawNucleotideLoci(), fragment.rawNucleotides()));

        fragment.addNucleotideInfo(9, "T", 30);
        assertTrue(fragment.containsNucleotide(9));
        assertTrue(fragment.validate());

        fragment.addNucleotideInfo(0, "A", 37);
        assertTrue(fragment.containsNucleotide(0));
        assertTrue(fragment.validate());
    }

    @Test
    public void testFragmentMerge()
    {
        String readId = "01";
        Read read = createReadRecord(readId);
        Fragment frag1 = new Fragment(
                read, GENE_A, Sets.newHashSet(GENE_A),
                Lists.newArrayList(1), Lists.newArrayList(30), Lists.newArrayList("A"));

        Fragment frag2 = new Fragment(
                read, GENE_B, Sets.newHashSet(GENE_B),
                Lists.newArrayList(1), Lists.newArrayList(30), Lists.newArrayList("A"));

        Fragment mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.genes().size());
        assertEquals(1, mergedFrag.nucleotideLoci().size());
        assertEquals(Integer.valueOf(1), mergedFrag.nucleotideLoci().get(0));
        assertEquals(1, mergedFrag.nucleotideQuality().size());
        assertEquals(1, mergedFrag.nucleotides().size());

        frag2 = new Fragment(
                read, GENE_A, Sets.newHashSet(GENE_A),
                Lists.newArrayList(0, 1, 2, 3),
                Lists.newArrayList(30, 30, 30, 30),
                Lists.newArrayList("A", "A", "A", "A"));

        mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.genes().size());
        assertEquals(4, mergedFrag.nucleotideLoci().size());
        assertEquals(Integer.valueOf(0), mergedFrag.nucleotideLoci().get(0));
        assertEquals(Integer.valueOf(1), mergedFrag.nucleotideLoci().get(1));
        assertEquals(4, mergedFrag.nucleotideQuality().size());
        assertEquals(4, mergedFrag.nucleotides().size());

        frag2 = new Fragment(
                read, GENE_C, Sets.newHashSet(GENE_C),
                Lists.newArrayList(3, 4, 5),
                Lists.newArrayList(30, 30, 30),
                Lists.newArrayList("A", "A", "A"));

        mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(3, mergedFrag.genes().size());
        assertEquals(6, mergedFrag.nucleotideLoci().size());
        assertEquals(Integer.valueOf(0), mergedFrag.nucleotideLoci().get(0));
        assertEquals(Integer.valueOf(3), mergedFrag.nucleotideLoci().get(3));
        assertEquals(Integer.valueOf(4), mergedFrag.nucleotideLoci().get(4));
        assertEquals(Integer.valueOf(5), mergedFrag.nucleotideLoci().get(5));
        assertEquals(6, mergedFrag.nucleotideQuality().size());
        assertEquals(6, mergedFrag.nucleotides().size());
    }

    private void assertRange(int expectedStart, int expectedEnd, List<Integer> victim)
    {
        assertEquals(expectedStart, victim.get(0).intValue());
        assertEquals(expectedEnd, victim.get(victim.size() - 1).intValue());
    }

    @Ignore
    @Test
    public void testLociLookup()
    {
        // test locus look-up for a fragment with typical number of nucleotides (~150 being a read within an exon, high value = 250)
        // relative search times: manual = 5, binary = 2, indexOf = 3
        int iterations = 1_000_000;

        int lociCount = 150;

        List<Integer> lociValues = Lists.newArrayListWithCapacity(lociCount);

        for(int i = 0; i < lociCount; ++i)
        {
            lociValues.add(i);
        }

        int[] searchLoci = {0, 10, 25, 50, 100, lociCount - 1};
        int searchValueIndex = 0;

        long startTime = System.nanoTime();

        for(int i = 0; i < iterations; ++i)
        {
            Fragment.findLociIndexManual(lociValues, searchLoci[searchValueIndex++]);

            if(searchValueIndex >= searchLoci.length)
                searchValueIndex = 0;
        }

        long sampleTime = System.nanoTime() - startTime;
        double sampleTimeMillis = sampleTime / NANOS_IN_MILLISECOND;

        LL_LOGGER.info(format("manual search time: %.6fms", sampleTimeMillis));

        searchValueIndex = 0;
        startTime = System.nanoTime();

        for(int i = 0; i < iterations; ++i)
        {
            Fragment.findLociInsertionIndex(lociValues, searchLoci[searchValueIndex++]);

            if(searchValueIndex >= searchLoci.length)
                searchValueIndex = 0;
        }

        sampleTime = System.nanoTime() - startTime;
        sampleTimeMillis = sampleTime / NANOS_IN_MILLISECOND;

        LL_LOGGER.info(format("binary search time: %.6fms", sampleTimeMillis));

        searchValueIndex = 0;
        startTime = System.nanoTime();

        for(int i = 0; i < iterations; ++i)
        {
            lociValues.indexOf(searchLoci[searchValueIndex++]);

            if(searchValueIndex >= searchLoci.length)
                searchValueIndex = 0;
        }

        sampleTime = System.nanoTime() - startTime;
        sampleTimeMillis = sampleTime / NANOS_IN_MILLISECOND;

        LL_LOGGER.info(format("indexOf search time: %.6fms", sampleTimeMillis));
    }
}
