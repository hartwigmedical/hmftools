package com.hartwig.hmftools.common.virusinterpreter;

import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VirusInterpreterTestFactory
{
    @NotNull
    public static ImmutableAnnotatedVirus.Builder annotatedVirusBuilder()
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
                .reported(false)
                .blacklisted(false)
                .virusDriverLikelihoodType(VirusLikelihoodType.HIGH);
    }
}
