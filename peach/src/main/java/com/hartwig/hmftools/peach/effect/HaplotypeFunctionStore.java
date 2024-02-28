package com.hartwig.hmftools.peach.effect;

import java.util.List;
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
        HaplotypeFunction matchingHaplotypeFunction = getMatchingHaplotypeFunctions(geneName, haplotypeName);
        return matchingHaplotypeFunction == null ? null : matchingHaplotypeFunction.function();
    }

    @Nullable
    private HaplotypeFunction getMatchingHaplotypeFunctions(@NotNull String geneName, @NotNull String haplotypeName)
    {
        List<HaplotypeFunction> matchingFunctions = haplotypeFunctions.stream()
                .filter(i -> i.geneName().equals(geneName) && i.haplotypeName().equals(haplotypeName))
                .collect(Collectors.toList());
        if(matchingFunctions.isEmpty())
        {
            return null;
        }
        else if(matchingFunctions.size() == 1)
        {
            return matchingFunctions.get(0);
        }
        else
        {
            String errorMessage = String.format("Multiple functions configured for haplotype %s for gene %s", haplotypeName, geneName);
            throw new IllegalStateException(errorMessage);
        }
    }
}
