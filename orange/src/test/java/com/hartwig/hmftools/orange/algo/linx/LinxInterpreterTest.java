package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.ImmutableLinxData;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LinxInterpreterTest {

    @Test
    public void canInterpretMinimalLinxData() {
        LinxInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(ImmutableLinxData.builder().build()));
    }

    @Test
    public void canFilterReportableGermlineBreakends() {
        LinxBreakend present = LinxTestFactory.breakendBuilder().junctionCopyNumber(1D).build();
        LinxBreakend absent = LinxTestFactory.breakendBuilder().junctionCopyNumber(0D).build();

        assertNull(LinxInterpreter.filterForPresenceInTumor(null));
        assertEquals(1, LinxInterpreter.filterForPresenceInTumor(Lists.newArrayList(present)).size());
        assertEquals(0, LinxInterpreter.filterForPresenceInTumor(Lists.newArrayList(absent)).size());
    }

    @Test
    public void canResolveHomozygousDisruptions() {
        LinxBreakend match = LinxTestFactory.breakendBuilder()
                .svId(1)
                .reportedDisruption(true)
                .junctionCopyNumber(2)
                .undisruptedCopyNumber(1)
                .type(StructuralVariantType.DUP)
                .build();

        LinxBreakend noPair = LinxTestFactory.breakendBuilder().svId(2).build();

        LinxBreakend notAHomPair = LinxTestFactory.breakendBuilder()
                .svId(1)
                .reportedDisruption(true)
                .junctionCopyNumber(1)
                .undisruptedCopyNumber(2)
                .type(StructuralVariantType.DUP)
                .build();

        LinxBreakend homPair = LinxTestFactory.breakendBuilder()
                .svId(1)
                .reportedDisruption(true)
                .junctionCopyNumber(2)
                .undisruptedCopyNumber(2)
                .type(StructuralVariantType.DUP)
                .build();

        assertNull(LinxInterpreter.resolveHomozygousDisruptions(null));
        assertTrue(LinxInterpreter.resolveHomozygousDisruptions(Lists.newArrayList(match, noPair)).isEmpty());
        assertTrue(LinxInterpreter.resolveHomozygousDisruptions(Lists.newArrayList(match, notAHomPair)).isEmpty());

        List<HomozygousDisruption> homozygousDisruptions = LinxInterpreter.resolveHomozygousDisruptions(Lists.newArrayList(match, homPair));
        assertEquals(1, homozygousDisruptions.size());
        HomozygousDisruption homozygousDisruption = homozygousDisruptions.get(0);
        assertEquals(match.chromosome(), homozygousDisruption.chromosome());
        assertEquals(match.chrBand(), homozygousDisruption.chromosomeBand());
        assertEquals(match.gene(), homozygousDisruption.gene());
        assertEquals(match.transcriptId(), homozygousDisruption.transcript());
        assertEquals(match.canonical(), homozygousDisruption.isCanonical());
    }

    @NotNull
    private static LinxInterpreter createTestInterpreter() {
        return new LinxInterpreter(Lists.newArrayList(), new KnownFusionCache());
    }
}