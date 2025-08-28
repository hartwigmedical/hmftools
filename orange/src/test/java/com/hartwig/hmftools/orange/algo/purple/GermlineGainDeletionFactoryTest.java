package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineGainDeletionFactoryTest
{
    private static final String TEST_GENE = "gene";
    private static final double EPSILON = 1.0E-2;

    @Test
    public void canFilterHetDeletion()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        GermlineDeletion reportableHet = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION);
        assertTrue(factory.getReportabilityMap(Lists.newArrayList(reportableHet), Lists.newArrayList()).isEmpty());
    }

    @Test
    public void canTransformReportableHomDeletionsToPartial()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        GermlineDeletion reportablePartialHom1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 400, 700);
        GermlineDeletion reportablePartialHom2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 800, 1000);
        List<GermlineDeletion> deletions = Lists.newArrayList(reportablePartialHom1, reportablePartialHom2);

        GeneCopyNumber partialLoss = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 1D, 4D);

        Map<PurpleGainDeletion, Boolean> map = factory.getReportabilityMap(deletions, Lists.newArrayList(partialLoss));
        PurpleGainDeletion gainDel = map.keySet().iterator().next();

        assertEquals(1, map.keySet().size());
        assertTrue(map.get(gainDel));
        assertEquals(CopyNumberInterpretation.PARTIAL_DEL, gainDel.interpretation());
        assertEquals(TEST_GENE, gainDel.gene());
        assertEquals(0, gainDel.minCopies(), EPSILON);
        assertEquals(4, gainDel.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformNonReportableHomDeletionToFull()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineDeletion reportableFullHom =
                GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HOM_DELETION, 0D, 100, 1200);
        GeneCopyNumber fullLoss = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 1D, 1D);

        Map<PurpleGainDeletion, Boolean> map = factory.getReportabilityMap(Lists.newArrayList(reportableFullHom), Lists.newArrayList(fullLoss));
        PurpleGainDeletion gainDel = map.keySet().iterator().next();

        assertEquals(1, map.keySet().size());
        assertFalse(map.get(gainDel));
        assertEquals(CopyNumberInterpretation.FULL_DEL, gainDel.interpretation());
        assertEquals(TEST_GENE, gainDel.gene());
        assertEquals(0, gainDel.minCopies(), EPSILON);
        assertEquals(0, gainDel.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformReportablePartialHomDeletionsToFullGeneLoss()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        GermlineDeletion reportablePartial1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 200, 300);
        GermlineDeletion reportablePartial2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 300, 500);
        GermlineDeletion reportablePartial3 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 500, 800);
        GermlineDeletion reportablePartial4 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.2D, 700, 2000);
        List<GermlineDeletion> deletions =
                Lists.newArrayList(reportablePartial1, reportablePartial2, reportablePartial3, reportablePartial4);

        GeneCopyNumber partialLoss = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 1D, 4D);

        Map<PurpleGainDeletion, Boolean> map = factory.getReportabilityMap(deletions, Lists.newArrayList(partialLoss));
        PurpleGainDeletion gainDel = map.keySet().iterator().next();

        assertEquals(1, map.keySet().size());
        assertTrue(map.get(gainDel));
        assertEquals(CopyNumberInterpretation.FULL_DEL, gainDel.interpretation());
        assertEquals(TEST_GENE, gainDel.gene());
        assertEquals(0, gainDel.minCopies(), EPSILON);
        assertEquals(0.2, gainDel.maxCopies(), EPSILON);
    }

    @NotNull
    private static GermlineGainDeletionFactory createTestFactory()
    {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        return new GermlineGainDeletionFactory(ensemblDataCache);
    }
}