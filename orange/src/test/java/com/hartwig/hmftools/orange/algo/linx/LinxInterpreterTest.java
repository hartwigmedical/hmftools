package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LinxInterpreterTest
{
    @Test
    public void canInterpretMinimalLinxData()
    {
        LinxInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(ImmutableLinxData.builder().build()));
    }

    @NotNull
    private static LinxInterpreter createTestInterpreter()
    {
        return new LinxInterpreter(TestEnsemblDataCacheFactory.loadTestCache());
    }
}