package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LinxInterpreterTest {

    @Test
    public void canFilterReportableGermlineBreakends() {
        LinxBreakend present = LinxTestFactory.breakendBuilder().junctionCopyNumber(1D).build();
        LinxBreakend absent = LinxTestFactory.breakendBuilder().junctionCopyNumber(0D).build();

        assertNull(LinxInterpreter.filterForPresenceInTumor(null));
        assertEquals(1, LinxInterpreter.filterForPresenceInTumor(Lists.newArrayList(present)).size());
        assertEquals(0, LinxInterpreter.filterForPresenceInTumor(Lists.newArrayList(absent)).size());
    }

    @Test
    public void canInterpretMinimalLinxData() {
        LinxInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(ImmutableLinxData.builder().build()));
    }

    @NotNull
    private static LinxInterpreter createTestInterpreter() {
        return new LinxInterpreter(Lists.newArrayList(), new KnownFusionCache());
    }
}