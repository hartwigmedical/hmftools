package com.hartwig.hmftools.orange.algo.virus;

import static java.lang.String.format;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.ImmutableVirus;
import com.hartwig.hmftools.datamodel.finding.Virus;
import com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.algo.util.FindingKeys;

import org.jetbrains.annotations.NotNull;

public final class VirusInterpreter
{
    @NotNull
    public static VirusInterpreterData interpret(@NotNull com.hartwig.hmftools.common.virus.VirusInterpreterData interpreterData)
    {
        List<VirusInterpreterEntry> allViruses = ConversionUtil.mapToList(
                filterBlacklistedViruses(interpreterData.allViruses()), VirusInterpreter::convert);

        List<Virus> driverViruses = allViruses.stream().filter(VirusInterpreterEntry::reported)
                .map(VirusInterpreter::toVirusFinding).toList();

        return com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterData.builder()
                .allViruses(allViruses)
                .driverViruses(driverViruses)
                .build();
    }

    @NotNull
    public static List<AnnotatedVirus> filterBlacklistedViruses(@NotNull List<AnnotatedVirus> allViruses)
    {
        Optional<AnnotatedVirus> virusWithBlacklistStatusUnknown =
                allViruses.stream().filter(virus -> virus.blacklisted() == null).findFirst();
        if(virusWithBlacklistStatusUnknown.isPresent())
        {
            throw new RuntimeException(format("Encountered virus '%s' with unknown blacklist status", virusWithBlacklistStatusUnknown.get()));
        }
        return allViruses.stream().filter(virus -> !virus.blacklisted()).collect(Collectors.toList());
    }

    @NotNull
    public static VirusInterpreterEntry convert(@NotNull com.hartwig.hmftools.common.virus.AnnotatedVirus annotatedVirus)
    {
        VirusType interpretation = annotatedVirus.interpretation();
        return ImmutableVirusInterpreterEntry.builder()
                .name(annotatedVirus.name())
                .qcStatus(VirusBreakendQCStatus.valueOf(annotatedVirus.qcStatus().name()))
                .integrations(annotatedVirus.integrations())
                .interpretation(interpretation != null ? VirusInterpretation.valueOf(interpretation.name()) : null)
                .percentageCovered(annotatedVirus.percentageCovered())
                .meanCoverage(annotatedVirus.meanCoverage())
                .expectedClonalCoverage(annotatedVirus.expectedClonalCoverage())
                .reported(annotatedVirus.reported())
                .driverLikelihood(VirusLikelihoodType.valueOf(annotatedVirus.virusDriverLikelihoodType().name()))
                .build();
    }

    @NotNull
    private static Virus toVirusFinding(@NotNull VirusInterpreterEntry virusInterpreterEntry)
    {
        if(!virusInterpreterEntry.reported())
        {
            throw new IllegalArgumentException("Cannot convert unreported virus to finding.");
        }
        return ImmutableVirus.builder()
                .findingKey(FindingKeys.virus(virusInterpreterEntry))
                .driverInterpretation(virusDriverInterpretation(virusInterpreterEntry.driverLikelihood()))
                .reportedStatus(ReportedStatus.REPORTED)
                .interpreterEntry(virusInterpreterEntry)
                .build();
    }

    @NotNull
    private static DriverInterpretation virusDriverInterpretation(@NotNull VirusLikelihoodType virusLikelihoodType)
    {
        return switch(virusLikelihoodType)
        {
            case LOW -> DriverInterpretation.LOW;
            case HIGH -> DriverInterpretation.HIGH;
            default -> throw new IllegalStateException("Unexpected virus likelihood type: " + virusLikelihoodType);
        };
    }
}
