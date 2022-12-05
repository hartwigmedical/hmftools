package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class DriverInterpretationTest {

    @Test
    public void canInterpretDriverLikelihood() {
        assertEquals(DriverInterpretation.HIGH, DriverInterpretation.interpret(0.9));
        assertEquals(DriverInterpretation.MEDIUM, DriverInterpretation.interpret(0.5));
        assertEquals(DriverInterpretation.LOW, DriverInterpretation.interpret(0.1));
    }

}