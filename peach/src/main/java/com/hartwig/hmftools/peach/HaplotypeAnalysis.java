package com.hartwig.hmftools.peach;

import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class HaplotypeAnalysis
{
    @NotNull
    private final Map<String, Integer> eventIdToCount;
    @NotNull
    private final List<HaplotypeCombination> haplotypeCombinations;

    public HaplotypeAnalysis(@NotNull Map<String, Integer> eventIdToCount, @NotNull List<HaplotypeCombination> haplotypeCombinations)
    {
        this.eventIdToCount = new HashMap<>(eventIdToCount);
        this.haplotypeCombinations = new ArrayList<>(haplotypeCombinations);
    }

    @Override
    public String toString()
    {
        return "HaplotypeAnalysis(" +
                "eventIdToCount=" + eventIdToCount +
                ", haplotypeCombinations=" + haplotypeCombinations +
                ')';
    }
}
