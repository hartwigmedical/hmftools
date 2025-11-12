package com.hartwig.hmftools.orange.report.finding;

import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.datamodel.finding.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.ImmutableVirus;
import com.hartwig.hmftools.datamodel.finding.Virus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;

public class VirusFactory {
    private VirusFactory() {
    }

    @NotNull
    public static List<Virus> convert(@NotNull Collection<VirusInterpreterEntry> virusInterpreterEntries)
    {
        return virusInterpreterEntries.stream().map(VirusFactory::convert)
                .toList();
    }

    private static Virus convert(VirusInterpreterEntry virusInterpreterEntry)
    {
        return ImmutableVirus.builder().from(virusInterpreterEntry)
                .findingKey(FindingKeys.findingKey(virusInterpreterEntry))
                .isReportable(virusInterpreterEntry.reported())
                .isCandidate(false)
                .driverInterpretation(convert(virusInterpreterEntry.driverLikelihood()))
                .build();
    }

    private static DriverInterpretation convert(VirusLikelihoodType virusLikelihoodType)
    {
        return switch (virusLikelihoodType)
        {
            case HIGH -> DriverInterpretation.HIGH;
            case LOW, UNKNOWN -> DriverInterpretation.LOW;
        };
    }
}
