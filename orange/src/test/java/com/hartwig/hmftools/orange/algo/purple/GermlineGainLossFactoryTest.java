package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineGainLossFactoryTest {

    private static final String TEST_GENE = "gene";

    @Test
    public void canFilterHetDeletion() {
        GermlineGainLossFactory factory = createTestFactory();

        GermlineDeletion reportableHet = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION);
        assertTrue(factory.mapDeletions(Lists.newArrayList(reportableHet), Lists.newArrayList()).isEmpty());
    }

//    @Test
//    public void canTransformHomDeletionToPartial() {
//        GermlineGainLossFactory factory = createTestFactory();
//
//        // Gene runs from 150 to 950
//        GermlineDeletion reportablePartialHom = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 400, 500);
//        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).maxCopyNumber(4D).build();
//
//        Map<PurpleGainLoss, GermlineDeletion> map = factory.mapDeletions(Lists.newArrayList(reportablePartialHom), Lists.newArrayList());
//        PurpleGainLoss gainLoss = map.keySet().iterator().next();
//    }

    @NotNull
    private static GermlineGainLossFactory createTestFactory() {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        return new GermlineGainLossFactory(ensemblDataCache);
    }
}