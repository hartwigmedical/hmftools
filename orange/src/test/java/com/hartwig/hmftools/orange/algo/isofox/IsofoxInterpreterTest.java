package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.orange.algo.linx.TestLinxFactory;

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
        return new IsofoxInterpreter(Maps.newHashMap(), TestLinxFactory.createMinimalTestLinxData());
    }
}