package com.hartwig.hmftools.common.virusinterpreter;

import com.hartwig.hmftools.common.viridian.VirusType;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendQCStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class VirusInterpreterTestFactory
{
    @NotNull
    public static VirusInterpreterData createEmptyData()
    {
        return ImmutableVirusInterpreterData.builder().build();
    }

    @NotNull
    public static VirusInterpreterData createMinimalData()
    {
        return createWithVirus(annotatedVirusBuilder().build());
    }

    @NotNull
    public static VirusInterpreterData createProperData()
    {
        return createWithVirus(annotatedVirusBuilder().interpretation(VirusType.MCV).expectedClonalCoverage(1.0).build());
    }

    @NotNull
    public static VirusInterpreterData createHHVInterpretationData()
    {
        return createWithVirus(annotatedVirusBuilder().interpretation(VirusType.HHV8).build());
    }

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

    @NotNull
    private static VirusInterpreterData createWithVirus(@NotNull AnnotatedVirus virus)
    {
        return ImmutableVirusInterpreterData.builder().addReportableViruses(virus).addAllViruses(virus).build();
    }
}
