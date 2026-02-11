package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.TestDataUtils.CYTO_BANDS;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.orange.TestDataUtils;

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
        return new LinxInterpreter(CYTO_BANDS);
    }
}