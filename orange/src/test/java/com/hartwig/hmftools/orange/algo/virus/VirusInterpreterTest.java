package com.hartwig.hmftools.orange.algo.virus;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.beust.jcommander.internal.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterTest
{
    @Test
    public void canFilterBlacklistedViruses()
    {
        List<AnnotatedVirus> filteredViruses =
                VirusInterpreter.filterOutBlacklistedViruses(Lists.newArrayList(createAnnotatedVirus("blacklisted virus", false, true), createAnnotatedVirus("non blacklisted virus", true, false)));
        assertEquals(1, filteredViruses.size());
        assertEquals(filteredViruses.iterator().next().name(), "non blacklisted virus");
    }

    @NotNull
    public static AnnotatedVirus createAnnotatedVirus(@NotNull String name, boolean reported, boolean blacklisted)
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

}
