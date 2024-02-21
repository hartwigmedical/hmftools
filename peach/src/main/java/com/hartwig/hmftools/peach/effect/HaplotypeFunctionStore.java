package com.hartwig.hmftools.peach.effect;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HaplotypeFunctionStore
{
    @NotNull
    private final List<HaplotypeFunction> haplotypeFunctions;

    public HaplotypeFunctionStore(@NotNull List<HaplotypeFunction> haplotypeFunctions)
    {
        this.haplotypeFunctions = haplotypeFunctions;
    }

    @Nullable
    public String getFunction(@NotNull String geneName, @NotNull String haplotypeName)
    {
        Set<String> configuredFunctions = haplotypeFunctions.stream()
                .filter(i -> i.geneName().equals(geneName) && i.haplotypeName().equals(haplotypeName))
                .map(HaplotypeFunction::function)
                .collect(Collectors.toSet());
        if(configuredFunctions.isEmpty())
        {
            return null;
        }
        else if(configuredFunctions.size() == 1)
        {
            return configuredFunctions.iterator().next();
        }
        else
        {
            String errorMessage =
                    String.format("Multiple different haplotype functions configured for haplotype %s of gene %s", haplotypeName, geneName);
            throw new IllegalStateException(errorMessage);
        }
    }
}
