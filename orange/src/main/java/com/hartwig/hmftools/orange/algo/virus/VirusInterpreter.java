package com.hartwig.hmftools.orange.algo.virus;

import static java.lang.String.format;

import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.OrangeConversion;

import org.jetbrains.annotations.NotNull;

public final class VirusInterpreter
{
    @NotNull
    public static VirusInterpreterData interpret(@NotNull com.hartwig.hmftools.common.virus.VirusInterpreterData interpreterData)
    {
        return com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterData.builder()
                .allViruses(ConversionUtil.mapToIterable(filterBlacklistedViruses(interpreterData.allViruses()), OrangeConversion::convert))
                .reportableViruses(ConversionUtil.mapToIterable(filterBlacklistedViruses(
                        interpreterData.reportableViruses()), OrangeConversion::convert))
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
}
