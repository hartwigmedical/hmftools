package com.hartwig.hmftools.lilac.fragment;

import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.LilacUtils.formRange;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.calcAminoAcidIndices;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.mergeFragments;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class FragmentsTest
{
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
    public void testFragmentMerge()
    {
        String readId = "01";
        Fragment frag1 = new Fragment(
                readId, "", GENE_A, Sets.newHashSet(GENE_A),
                Lists.newArrayList(1), Lists.newArrayList(30), Lists.newArrayList("A"));

        Fragment frag2 = new Fragment(
                readId, "", GENE_B, Sets.newHashSet(GENE_B),
                Lists.newArrayList(1), Lists.newArrayList(30), Lists.newArrayList("A"));

        Fragment mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.getGenes().size());
        assertEquals(1, mergedFrag.getNucleotideLoci().size());
        assertEquals(Integer.valueOf(1), mergedFrag.getNucleotideLoci().get(0));
        assertEquals(1, mergedFrag.getNucleotideQuality().size());
        assertEquals(1, mergedFrag.getNucleotides().size());

        frag2 = new Fragment(
                readId, "", GENE_A, Sets.newHashSet(GENE_A),
                Lists.newArrayList(0, 1, 2, 3),
                Lists.newArrayList(30, 30, 30, 30),
                Lists.newArrayList("A", "A", "A", "A"));

        mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(2, mergedFrag.getGenes().size());
        assertEquals(4, mergedFrag.getNucleotideLoci().size());
        assertEquals(Integer.valueOf(0), mergedFrag.getNucleotideLoci().get(0));
        assertEquals(Integer.valueOf(1), mergedFrag.getNucleotideLoci().get(1));
        assertEquals(4, mergedFrag.getNucleotideQuality().size());
        assertEquals(4, mergedFrag.getNucleotides().size());

        frag2 = new Fragment(
                readId, "", GENE_C, Sets.newHashSet(GENE_C),
                Lists.newArrayList(3, 4, 5),
                Lists.newArrayList(30, 30, 30),
                Lists.newArrayList("A", "A", "A"));

        mergedFrag = mergeFragments(frag1, frag2);
        assertTrue(frag1.validate());
        assertEquals(3, mergedFrag.getGenes().size());
        assertEquals(6, mergedFrag.getNucleotideLoci().size());
        assertEquals(Integer.valueOf(0), mergedFrag.getNucleotideLoci().get(0));
        assertEquals(Integer.valueOf(3), mergedFrag.getNucleotideLoci().get(3));
        assertEquals(Integer.valueOf(4), mergedFrag.getNucleotideLoci().get(4));
        assertEquals(Integer.valueOf(5), mergedFrag.getNucleotideLoci().get(5));
        assertEquals(6, mergedFrag.getNucleotideQuality().size());
        assertEquals(6, mergedFrag.getNucleotides().size());
    }

    private void assertRange(int expectedStart, int expectedEnd, List<Integer> victim)
    {
        assertEquals(expectedStart, victim.get(0).intValue());
        assertEquals(expectedEnd, victim.get(victim.size() - 1).intValue());
    }

}
