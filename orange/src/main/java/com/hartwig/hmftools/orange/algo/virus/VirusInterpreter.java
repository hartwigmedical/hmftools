package com.hartwig.hmftools.orange.algo.virus;

import java.util.List;
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
        return allViruses.stream().filter(virus -> !virus.blacklisted()).collect(Collectors.toList());
    }
}
