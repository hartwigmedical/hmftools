package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;

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
        return new PurpleInterpreter(Lists.newArrayList(), ChordTestFactory.createMinimalTestChordAnalysis());
    }
}