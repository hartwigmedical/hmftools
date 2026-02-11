package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.orange.TestDataUtils;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineAmpDelUtilTest
{
    private static final String TEST_GENE = "gene";
    private static final String TEST_GENE2 = "gene2";  // not compatible with test Ensembl data cache
    private static final String TEST_GENE3 = "gene3";  // not compatible with test Ensembl data cache
    private static final double EPSILON = 1.0E-2;

    @Test
    public void canCalculateMaxCopyNumberForPartialDeletionsHetInTumorWithHighSomaticCopyNumber()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel partialDeletionHetInTumor1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.1D, 300, 400);
        GermlineAmpDel partialDeletionHetInTumor2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 0.9D, 600, 800);
        List<GermlineAmpDel> deletions = Lists.newArrayList(partialDeletionHetInTumor1, partialDeletionHetInTumor2);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 2D);
        double maxCopyNumber = GermlineGainDeletionUtil.getSomaticMaxCopyNumber(deletions, geneCopyNumber, transcript);
        assertEquals(2.0, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForPartialDeletionsHetInTumorWithLowSomaticCopyNumber()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel partialDeletionHetInTumor1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.1D, 300, 400);
        GermlineAmpDel partialDeletionHetInTumor2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 0.9D, 600, 800);
        List<GermlineAmpDel> deletions = Lists.newArrayList(partialDeletionHetInTumor1, partialDeletionHetInTumor2);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 0.8D);
        double maxCopyNumber = GermlineGainDeletionUtil.getSomaticMaxCopyNumber(deletions, geneCopyNumber, transcript);
        assertEquals(1.1, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForPartialDeletionsHomInTumor()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel partialDeletionHomInTumor1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 300, 400);
        GermlineAmpDel partialDeletionHomInTumor2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.2D, 600, 800);
        List<GermlineAmpDel> deletions = Lists.newArrayList(partialDeletionHomInTumor1, partialDeletionHomInTumor2);

        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 2D);
        double maxCopyNumber = GermlineGainDeletionUtil.getSomaticMaxCopyNumber(deletions, geneCopyNumber, transcript);
        assertEquals(2.0, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForFullDeletionHetInTumor()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel fullDeletionHetInTumor =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.1D, 100, 1000);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 2D);
        double maxCopyNumber =
                GermlineGainDeletionUtil.getSomaticMaxCopyNumber(Lists.newArrayList(fullDeletionHetInTumor), geneCopyNumber, transcript);
        assertEquals(1.1, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForFullDeletionHomInTumor()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel fullDeletionHomInTumor =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 2D);
        double maxCopyNumber =
                GermlineGainDeletionUtil.getSomaticMaxCopyNumber(Lists.newArrayList(fullDeletionHomInTumor), geneCopyNumber, transcript);
        assertEquals(0.1, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForFullDeletionWithNegativeCNHomInTumor()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel fullDeletionHomInTumor =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, -0.1D, 100, 1000);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 2D);
        double maxCopyNumber =
                GermlineGainDeletionUtil.getSomaticMaxCopyNumber(Lists.newArrayList(fullDeletionHomInTumor), geneCopyNumber, transcript);
        assertEquals(0, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForPartialDeletionsTogetherFullDeletion()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();
        GermlineAmpDel partialDeletionHetInTumor1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.1D, 200, 700);
        GermlineAmpDel partialDeletionHomInTumor2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.2D, 701, 1000);
        List<GermlineAmpDel> deletions = Lists.newArrayList(partialDeletionHetInTumor1, partialDeletionHomInTumor2);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 0, 2D);
        double maxCopyNumber = GermlineGainDeletionUtil.getSomaticMaxCopyNumber(deletions, geneCopyNumber, transcript);
        assertEquals(1.1, maxCopyNumber, EPSILON);
    }

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
    public void testDeletionCoversTranscriptSimple()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();

        GermlineAmpDel fullDeletion = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        GermlineAmpDel firstHalfDeletion =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.2D, 100, 500);
        GermlineAmpDel secondHalfDeletion =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.3D, 501, 1000);

        assertTrue(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(fullDeletion), transcript));
        assertFalse(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(firstHalfDeletion), transcript));
        assertFalse(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(secondHalfDeletion), transcript));
        assertTrue(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(fullDeletion, firstHalfDeletion), transcript));
        assertTrue(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(firstHalfDeletion, secondHalfDeletion), transcript));
        assertTrue(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(fullDeletion, firstHalfDeletion, secondHalfDeletion), transcript));
    }

    @Test
    public void testDeletionCoversTranscriptDifficult()
    {
        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        TranscriptData transcript = getTestTranscript();

        GermlineAmpDel deletion1 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 150);
        GermlineAmpDel deletion2 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.2D, 200, 300);
        GermlineAmpDel deletion3 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.D, 275, 325);
        GermlineAmpDel deletion4 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.D, 302, 700);
        GermlineAmpDel deletion5 = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 3.1D, 701, 1000);

        assertTrue(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(deletion1, deletion2, deletion3, deletion4, deletion5), transcript));
        assertFalse(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(deletion1, deletion2, deletion4, deletion5), transcript));
        assertFalse(GermlineGainDeletionUtil.deletionsCoverTranscript(Lists.newArrayList(deletion1), transcript));
    }

    @Test
    public void testGetCanonicalTranscript()
    {
        EnsemblDataCache ensemblDataCache = TestDataUtils.loadTestCache();
        TranscriptData canonicalTranscript = GermlineGainDeletionUtil.findCanonicalTranscript(TEST_GENE, ensemblDataCache);

        assertEquals("gene_id", canonicalTranscript.GeneId);
        assertTrue(canonicalTranscript.IsCanonical);
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

    @NotNull
    private static TranscriptData getTestTranscript()
    {
        EnsemblDataCache ensemblDataCache = TestDataUtils.loadTestCache();
        return ensemblDataCache.getCanonicalTranscriptData("gene_id");
    }
}