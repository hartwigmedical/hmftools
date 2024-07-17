package com.hartwig.hmftools.orange.algo.virus;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.beust.jcommander.internal.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.common.virus.VirusTestFactory;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterTest
{
    private static final double EPSILON = 0.01;

    @Test
    public void shouldConvertAnnotatedVirus()
    {
        assertTrue(VirusInterpreter.interpret(VirusTestFactory.createEmptyData()).allViruses().isEmpty());
        assertEqualsValue(VirusTestFactory.createMinimalData(), VirusInterpreter.interpret(VirusTestFactory.createMinimalData()));
        assertEqualsValue(VirusTestFactory.createProperData(), VirusInterpreter.interpret(VirusTestFactory.createProperData()));
        assertEqualsValue(VirusTestFactory.createHHVInterpretationData(), VirusInterpreter.interpret(VirusTestFactory.createHHVInterpretationData()));
    }

    @Test
    public void canFilterBlacklistedViruses()
    {
        List<AnnotatedVirus> filteredViruses =
                VirusInterpreter.filterBlacklistedViruses(Lists.newArrayList(createAnnotatedVirus("blacklisted virus", false, true), createAnnotatedVirus("non blacklisted virus", true, false)));
        assertEquals(1, filteredViruses.size());
        assertEquals(filteredViruses.iterator().next().name(), "non blacklisted virus");
    }

    @NotNull
    private static AnnotatedVirus createAnnotatedVirus(@NotNull String name, boolean reported, boolean blacklisted)
    {
        return ImmutableAnnotatedVirus.builder()
                .taxid(0)
                .name(name)
                .interpretation(null)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .expectedClonalCoverage(1.0)
                .reported(reported)
                .blacklisted(blacklisted)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH).build();
    }

    private static void assertEqualsValue(@NotNull com.hartwig.hmftools.common.virus.VirusInterpreterData input,
            @NotNull VirusInterpreterData converted)
    {
        assertEqualsValue(input.allViruses(), converted.allViruses());
        assertEqualsValue(input.reportableViruses(), converted.reportableViruses());
    }

    private static void assertEqualsValue(@NotNull List<com.hartwig.hmftools.common.virus.AnnotatedVirus> input,
            @NotNull List<VirusInterpreterEntry> converted)
    {
        assertEquals(converted.size(), input.size());
        for(int i = 0; i < input.size(); i++)
        {
            assertEqualsValue(input.get(i), converted.get(i));
        }
    }

    private static void assertEqualsValue(@NotNull com.hartwig.hmftools.common.virus.AnnotatedVirus input,
            @NotNull VirusInterpreterEntry converted)
    {
        assertEquals(input.name(), converted.name());
        assertEquals(input.qcStatus().name(), converted.qcStatus().name());
        assertEquals(input.integrations(), converted.integrations());
        assertEquals(input.percentageCovered(), converted.percentageCovered(), EPSILON);
        assertEquals(input.meanCoverage(), converted.meanCoverage(), EPSILON);
        assertEquals(input.expectedClonalCoverage(), converted.expectedClonalCoverage());
        assertEquals(input.reported(), converted.reported());
        assertEquals(input.virusDriverLikelihoodType().name(), converted.driverLikelihood().name());
    }
}
