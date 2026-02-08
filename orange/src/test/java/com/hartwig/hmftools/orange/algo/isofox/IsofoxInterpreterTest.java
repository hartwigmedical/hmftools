package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class IsofoxInterpreterTest
{
    @Test
    public void canInterpretMinimalIsofoxData()
    {
        IsofoxInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(IsofoxTestFactory.createMinimalIsofoxTestData()));
    }

    @NotNull
    private static IsofoxInterpreter createTestInterpreter()
    {
        return new IsofoxInterpreter(Lists.newArrayList(), TestLinxInterpretationFactory.createMinimalTestLinxData());
    }
}