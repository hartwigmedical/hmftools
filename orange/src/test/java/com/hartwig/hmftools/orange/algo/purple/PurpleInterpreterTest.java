package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutablePurpleData;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleInterpreterTest {

    private static final String TEST_GENE = "gene";

    @Test
    public void canCreateGermlineGainsLosses() {
        // Gene is needed to be able to match with ensembl test data
        GermlineDeletion hetReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1);
        GermlineDeletion homReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 2);
        GermlineDeletion homUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HOM_DELETION, 3);

        PurpleData purple = ImmutablePurpleData.builder()
                .from(PurpleTestFactory.createMinimalTestPurpleData())
                .addAllSomaticGeneCopyNumbers(GeneCopyNumberTestFactory.builder().chromosome("1").geneName(TEST_GENE).build())
                .addAllGermlineDeletions(hetReported, homReported, homUnreported)
                .addReportableGermlineDeletions(hetReported, homReported)
                .build();

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleInterpretedData interpreted = interpreter.interpret(purple);

        assertEquals(2, interpreted.allGermlineGainsLosses().size());
        assertEquals(1, interpreted.reportableGermlineGainsLosses().size());
    }

    @Test
    public void canInterpretMinimalPurpleData() {
        PurpleInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(PurpleTestFactory.createMinimalTestPurpleData()));
    }

    @NotNull
    private static PurpleInterpreter createTestInterpreter() {
        return createInterpreter(TestEnsemblDataCacheFactory.createDummyCache());
    }

    @NotNull
    private static PurpleInterpreter createRealInterpreter() {
        return createInterpreter(TestEnsemblDataCacheFactory.loadTestCache());
    }

    private static PurpleInterpreter createInterpreter(@NotNull EnsemblDataCache ensemblDataCache) {
        PurpleVariantFactory purpleVariantFactory = new PurpleVariantFactory(new PaveAlgo(ensemblDataCache));
        GermlineGainLossFactory germlineGainLossFactory = new GermlineGainLossFactory(ensemblDataCache);
        return new PurpleInterpreter(purpleVariantFactory, germlineGainLossFactory, Lists.newArrayList(), null);
    }
}