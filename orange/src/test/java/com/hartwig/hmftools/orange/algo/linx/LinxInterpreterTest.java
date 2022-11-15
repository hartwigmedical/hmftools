package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class LinxInterpreterTest {

    @Test
    public void canInterpretMinimalLinxData() {
        assertNotNull(LinxInterpreter.interpret(ImmutableLinxData.builder().build(),
                Lists.newArrayList(),
                new KnownFusionCache()));
    }
}