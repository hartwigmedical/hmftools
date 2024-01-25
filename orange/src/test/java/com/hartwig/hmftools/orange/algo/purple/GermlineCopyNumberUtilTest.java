package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineCopyNumberUtilTest
{
    private static final String TEST_GENE = "gene";
    private static final double EPSILON = 1.0E-2;

    @Test
    public void canCalculateMaxCopyNumberForPartialDeletionHetInTumor()
    {
        // Gene runs from 150 to 950
        GermlineDeletion partialDeletionHetInTumor = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HET_DELETION, 1.1D, 300, 400);
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).maxCopyNumber(2D).build();
        double maxCopyNumber =
                GermlineCopyNumberUtil.getSomaticMaxCopyNumber(partialDeletionHetInTumor, Lists.newArrayList(geneCopyNumber), ensemblDataCache);
        assertEquals(2.0, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForPartialDeletionHomInTumor()
    {
        // Gene runs from 150 to 950
        GermlineDeletion partialDeletionHomInTumor = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 300, 400);
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).maxCopyNumber(2D).build();
        double maxCopyNumber =
                GermlineCopyNumberUtil.getSomaticMaxCopyNumber(partialDeletionHomInTumor, Lists.newArrayList(geneCopyNumber), ensemblDataCache);
        assertEquals(2.0, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForFullDeletionHetInTumor()
    {
        // Gene runs from 150 to 950
        GermlineDeletion fullDeletionHetInTumor = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HET_DELETION, 1.1D, 100, 1000);
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).maxCopyNumber(2D).build();
        double maxCopyNumber =
                GermlineCopyNumberUtil.getSomaticMaxCopyNumber(fullDeletionHetInTumor, Lists.newArrayList(geneCopyNumber), ensemblDataCache);
        assertEquals(1.1, maxCopyNumber, EPSILON);
    }

    @Test
    public void canCalculateMaxCopyNumberForFullDeletionHomInTumor()
    {
        // Gene runs from 150 to 950
        GermlineDeletion fullDeletionHomInTumor = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).maxCopyNumber(2D).build();
        double maxCopyNumber =
                GermlineCopyNumberUtil.getSomaticMaxCopyNumber(fullDeletionHomInTumor, Lists.newArrayList(geneCopyNumber), ensemblDataCache);
        assertEquals(0.1, maxCopyNumber, EPSILON);
    }

    @Test
    public void testGetSomaticMinCopyNumber()
    {
        GermlineDeletion deletion = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        assertEquals(0.1, GermlineCopyNumberUtil.getSomaticMinCopyNumber(deletion), EPSILON);
    }

    @Test
    public void testGetSomaticMinCopyNumberWithNegativeCopyNumber()
    {
        GermlineDeletion deletionWithNegativeCopyNumber = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, -0.3D, 100, 1000);
        assertEquals(0., GermlineCopyNumberUtil.getSomaticMinCopyNumber(deletionWithNegativeCopyNumber), EPSILON);
    }

    @Test
    public void testDeletionCoversTranscript()
    {
        // Gene runs from 150 to 950, coding from 300 to 800
        TranscriptData transcript = getTestTranscript();

        GermlineDeletion fullDeletion = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 1000);
        assertTrue(GermlineCopyNumberUtil.deletionCoversTranscript(fullDeletion, transcript));

        GermlineDeletion partialStartDeletion = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 100, 500);
        assertFalse(GermlineCopyNumberUtil.deletionCoversTranscript(partialStartDeletion, transcript));

        GermlineDeletion partialEndDeletion = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 400, 1000);
        assertFalse(GermlineCopyNumberUtil.deletionCoversTranscript(partialEndDeletion, transcript));

        GermlineDeletion partialInnerDeletion = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 400, 500);
        assertFalse(GermlineCopyNumberUtil.deletionCoversTranscript(partialInnerDeletion, transcript));

        GermlineDeletion codingDeletion = GermlineDeletionTestFactory.create(
                TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 200, 925);
        assertFalse(GermlineCopyNumberUtil.deletionCoversTranscript(codingDeletion, transcript));
    }

    @Test
    public void testGetCanonicalTranscript()
    {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        TranscriptData canonicalTranscript = GermlineCopyNumberUtil.findCanonicalTranscript(TEST_GENE, ensemblDataCache);

        assertEquals("gene_id", canonicalTranscript.GeneId);
        assertTrue(canonicalTranscript.IsCanonical);
    }

    @NotNull
    private static TranscriptData getTestTranscript()
    {
        return new TranscriptData(
                1,
                "trans",
                "gene_id",
                true,
                Strand.NEG_STRAND,
                150,
                950,
                300,
                800,
                "protein_coding"
        );
    }
}