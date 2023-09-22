package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineGainLossFactoryTest
{
    private static final String TEST_GENE = "gene";
    private static final double EPSILON = 1.0E-2;

    @Test
    public void canFilterHetDeletion()
    {
        GermlineGainLossFactory factory = createTestFactory();

        GermlineDeletion reportableHet = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION);
        assertTrue(factory.mapDeletions(Lists.newArrayList(reportableHet), Lists.newArrayList()).isEmpty());
    }

    @Test
    public void canTransformHomDeletionToPartial()
    {
        GermlineGainLossFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineDeletion reportablePartialHom =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 400, 500);
        GeneCopyNumber partialLoss = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(0D).maxCopyNumber(4D).build();

        Map<PurpleGainLoss, GermlineDeletion> map =
                factory.mapDeletions(Lists.newArrayList(reportablePartialHom), Lists.newArrayList(partialLoss));
        PurpleGainLoss gainLoss = map.keySet().iterator().next();

        assertEquals(reportablePartialHom, map.get(gainLoss));
        assertEquals(CopyNumberInterpretation.PARTIAL_LOSS, gainLoss.interpretation());
        assertEquals(TEST_GENE, gainLoss.gene());
        assertEquals(0, gainLoss.minCopies(), EPSILON);
        assertEquals(4, gainLoss.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformHomDeletionToFull()
    {
        GermlineGainLossFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineDeletion reportableFullHom =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 100, 1200);
        GeneCopyNumber fullLoss = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(0D).maxCopyNumber(0D).build();

        Map<PurpleGainLoss, GermlineDeletion> map =
                factory.mapDeletions(Lists.newArrayList(reportableFullHom), Lists.newArrayList(fullLoss));
        PurpleGainLoss gainLoss = map.keySet().iterator().next();

        assertEquals(reportableFullHom, map.get(gainLoss));
        assertEquals(CopyNumberInterpretation.FULL_LOSS, gainLoss.interpretation());
        assertEquals(TEST_GENE, gainLoss.gene());
        assertEquals(0, gainLoss.minCopies(), EPSILON);
        assertEquals(0, gainLoss.maxCopies(), EPSILON);
    }

    @NotNull
    private static GermlineGainLossFactory createTestFactory()
    {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        return new GermlineGainLossFactory(ensemblDataCache);
    }
}