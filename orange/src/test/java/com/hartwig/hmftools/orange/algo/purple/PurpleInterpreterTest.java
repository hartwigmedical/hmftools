package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleInterpreterTest {

    @Test
    public void canInterpretMinimalPurpleData() {
        PurpleInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(PurpleTestFactory.createMinimalTestPurpleData()));
    }

    @NotNull
    private static PurpleInterpreter createTestInterpreter() {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.createDummyCache();
        PurpleVariantFactory purpleVariantFactory = new PurpleVariantFactory(new PaveAlgo(ensemblDataCache));
        GermlineGainLossFactory germlineGainLossFactory = new GermlineGainLossFactory(ensemblDataCache);
        return new PurpleInterpreter(purpleVariantFactory, germlineGainLossFactory, Lists.newArrayList(), null);
    }
}