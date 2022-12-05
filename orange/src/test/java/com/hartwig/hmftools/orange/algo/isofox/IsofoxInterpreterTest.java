package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertNotNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.isofox.IsofoxTestFactory;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;

import org.junit.Test;

public class IsofoxInterpreterTest {

    @Test
    public void canInterpretMinimalIsofoxData() {
        assertNotNull(IsofoxInterpreter.interpret(IsofoxTestFactory.createMinimalIsofoxTestData(),
                TestLinxInterpretationFactory.createMinimalTestLinxData(),
                Lists.newArrayList(),
                new KnownFusionCache()));
    }
}