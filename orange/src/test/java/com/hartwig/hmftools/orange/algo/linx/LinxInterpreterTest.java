package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;

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
        return new LinxInterpreter(Lists.newArrayList(), new KnownFusionCache());
    }
}