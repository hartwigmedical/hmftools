package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.junit.Test;

public class GermlineAmpDelUtilTest
{
    private static final String TEST_GENE = "gene";
    private static final String TEST_GENE2 = "gene2";  // not compatible with test Ensembl data cache
    private static final String TEST_GENE3 = "gene3";  // not compatible with test Ensembl data cache
    private static final double EPSILON = 1.0E-2;

    @Test
    public void testGetSomaticMinCopyNumber()
    {
        GermlineAmpDel deletion1 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.2D, 100, 1000);
        GermlineAmpDel deletion2 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        assertEquals(0.1, GermlineGainDeletionUtil.getSomaticMinCopyNumber(Lists.newArrayList(deletion1, deletion2)), EPSILON);
    }

    @Test
    public void testGetSomaticMinCopyNumberWithNegativeCopyNumber()
    {
        GermlineAmpDel deletion1 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, -0.3D, 100, 1000);
        GermlineAmpDel deletion2 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        assertEquals(0, GermlineGainDeletionUtil.getSomaticMinCopyNumber(Lists.newArrayList(deletion1, deletion2)), EPSILON);
    }

    @Test
    public void testGetChromosome()
    {
        GermlineAmpDel deletion1 = GermlineDeletionTestFactory.create(TEST_GENE, true, "1", "q20.30");
        GermlineAmpDel deletion2 = GermlineDeletionTestFactory.create(TEST_GENE, true, "1", "q20.30");
        GermlineAmpDel deletion3 = GermlineDeletionTestFactory.create(TEST_GENE, true, "1", "q20.30");
        List<GermlineAmpDel> deletions = Lists.newArrayList(deletion1, deletion2, deletion3);
        assertEquals("1", GermlineGainDeletionUtil.getChromosome(deletions));
    }

    @Test
    public void testGetChromosomeBand()
    {
        GermlineAmpDel deletion1 = GermlineDeletionTestFactory.create(TEST_GENE, true, "1", "q20.30");
        GermlineAmpDel deletion2 = GermlineDeletionTestFactory.create(TEST_GENE, true, "1", "q20.30");
        List<GermlineAmpDel> deletions = Lists.newArrayList(deletion1, deletion2);
        assertEquals("q20.30", GermlineGainDeletionUtil.getChromosomeBand(deletions));
    }

    @Test
    public void testFindCopyNumberForGene()
    {
        GeneCopyNumber geneCopyNumberToFind = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE);
        GeneCopyNumber otherGeneCopyNumber1 = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE2);
        GeneCopyNumber otherGeneCopyNumber2 = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE3);
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList(otherGeneCopyNumber1, geneCopyNumberToFind, otherGeneCopyNumber2);
        GeneCopyNumber foundGeneCopyNumber = GermlineGainDeletionUtil.findGeneCopyNumberForGene(TEST_GENE, geneCopyNumbers);
        assertEquals(geneCopyNumberToFind, foundGeneCopyNumber);
    }
}