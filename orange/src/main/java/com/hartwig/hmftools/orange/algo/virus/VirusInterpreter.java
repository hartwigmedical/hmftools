package com.hartwig.hmftools.orange.algo.virus;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.virus.AnnotatedVirus;

import org.jetbrains.annotations.NotNull;

public class VirusInterpreter
{
    @NotNull
    public static List<AnnotatedVirus> filterBlacklistedViruses(@NotNull List<AnnotatedVirus> allViruses)
    {
        return allViruses.stream().filter(virus -> !virus.blacklisted()).collect(Collectors.toList());
    }
}
