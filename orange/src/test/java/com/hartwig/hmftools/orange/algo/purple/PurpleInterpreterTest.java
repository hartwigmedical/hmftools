package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordTestFactory;
import com.hartwig.hmftools.common.purple.PurpleTestFactory;

import org.junit.Test;

public class PurpleInterpreterTest {

    @Test
    public void canInterpretMinimalPurpleData() {
        assertNotNull(PurpleInterpreter.interpret(PurpleTestFactory.createMinimalTestPurpleData(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                ChordTestFactory.createMinimalTestChordAnalysis()));
    }
}