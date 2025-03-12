package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.mergeFragments;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.buildSamRecord;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createReadRecord;
import static com.hartwig.hmftools.lilac.read.ReadRecord.create;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.lilac.read.ReadRecord;

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

        ReadRecord readRecord = create(codingRegion, record, true, true);
        assertNotNull(readRecord);

        assertEquals(190, readRecord.PositionStart);
        assertEquals(269, readRecord.PositionEnd);
        assertEquals(10, readRecord.SoftClippedStart);
        assertEquals(10, readRecord.SoftClippedEnd);

        // now restricted at the lower 3' end
        record.setInferredInsertSize(70);
        record.setReadNegativeStrandFlag(true);

        readRecord = create(codingRegion, record, true, true);
        assertEquals(200, readRecord.PositionStart);
        assertEquals(269, readRecord.PositionEnd);
        assertEquals(0, readRecord.SoftClippedStart);
        assertEquals(10, readRecord.SoftClippedEnd);

        record.setReadNegativeStrandFlag(false);

        readRecord = create(codingRegion, record, true, true);
        assertEquals(190, readRecord.PositionStart);
        assertEquals(259, readRecord.PositionEnd);
        assertEquals(10, readRecord.SoftClippedStart);
        assertEquals(0, readRecord.SoftClippedEnd);
    }

    @Test
    public void testFragmentAddNucleotide()
    {
        List<Integer> indices = Lists.newArrayList(1, 2, 3, 6, 7, 8);
        List<Integer> qualities = Lists.newArrayList(37, 25, 37, 37, 37, 25);
        List<String> nucleotides = Lists.newArrayList("A", "G", "T", "C", "A", "G");

        Fragment fragment = new Fragment(createReadRecord("01"), HLA_A, Sets.newHashSet(HLA_A), indices, qualities, nucleotides);

        fragment.qualityFilter(30);

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
        ReadRecord readRecord = createReadRecord(readId);
        Fragment frag1 = new Fragment(
                readRecord, GENE_A, Sets.newHashSet(GENE_A),
                Lists.newArrayList(1), Lists.newArrayList(30), Lists.newArrayList("A"));

        Fragment frag2 = new Fragment(
                readRecord, GENE_B, Sets.newHashSet(GENE_B),
                Lists.newArrayList(1), Lists.newArrayList(30), Lists.newArrayList("A"));

        Fragment mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.genes().size());
        assertEquals(1, mergedFrag.nucleotideLoci().size());
        assertEquals(Integer.valueOf(1), mergedFrag.nucleotideLoci().get(0));
        assertEquals(1, mergedFrag.nucleotideQuality().size());
        assertEquals(1, mergedFrag.nucleotides().size());

        frag2 = new Fragment(
                readRecord, GENE_A, Sets.newHashSet(GENE_A),
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
                readRecord, GENE_C, Sets.newHashSet(GENE_C),
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
}
