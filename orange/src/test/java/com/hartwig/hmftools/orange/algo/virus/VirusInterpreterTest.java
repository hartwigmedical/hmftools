package com.hartwig.hmftools.orange.algo.virus;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.beust.jcommander.internal.Lists;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VirusInterpreterTest
{
    @Test
    public void canFilterBlacklistedViruses()
    {
        List<AnnotatedVirus> filteredViruses =
                VirusInterpreter.filterOutBlacklistedViruses(Lists.newArrayList(createAnnotatedVirus(true), createAnnotatedVirus(false)));
        assertEquals(1, filteredViruses.size());
    }

    @NotNull
    public static AnnotatedVirus createAnnotatedVirus(boolean blacklisted)
    {
        return ImmutableAnnotatedVirus.builder()
                .taxid(0)
                .name(Strings.EMPTY)
                .interpretation(null)
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .percentageCovered(1.0)
                .meanCoverage(1.0)
                .expectedClonalCoverage(1.0)
                .reported(true)
                .blacklisted(blacklisted)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH).build();
    }

}
