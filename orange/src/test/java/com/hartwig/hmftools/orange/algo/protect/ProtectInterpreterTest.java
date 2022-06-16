package com.hartwig.hmftools.orange.algo.protect;

import static org.junit.Assert.assertNotNull;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class ProtectInterpreterTest {

    @Test
    public void canInterpretMinimalProtectData() {
        assertNotNull(ProtectInterpreter.interpret(Lists.newArrayList()));
    }
}