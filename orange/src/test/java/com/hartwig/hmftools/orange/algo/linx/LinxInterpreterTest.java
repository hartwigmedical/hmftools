package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;
import com.hartwig.hmftools.orange.algo.purple.PurpleData;
import com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LinxInterpreterTest
{
    @Test
    public void canInterpretMinimalLinxData()
    {
        PurpleData purpleData = PurpleTestFactory.createMinimalTestPurpleData();
        LinxInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(ImmutableLinxData.builder().build(), purpleData.purityContext().qc()));
    }

    @NotNull
    private static LinxInterpreter createTestInterpreter()
    {
        return new LinxInterpreter(TestEnsemblDataCacheFactory.loadTestCache());
    }
}