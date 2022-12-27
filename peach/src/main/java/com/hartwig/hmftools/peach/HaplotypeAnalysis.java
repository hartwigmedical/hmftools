package com.hartwig.hmftools.peach;

import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

public class HaplotypeAnalysis
{
    @NotNull
    private final Map<String, Integer> eventIdToCount;
    @NotNull
    private final List<HaplotypeCombination> haplotypeCombinations;
    @NotNull
    private final String wildTypeHaplotypeName;

    public HaplotypeAnalysis(
            @NotNull Map<String, Integer> eventIdToCount,
            @NotNull List<HaplotypeCombination> haplotypeCombinations,
            @NotNull String wildTypeHaplotypeName
    )
    {
        this.eventIdToCount = new HashMap<>(eventIdToCount);
        this.haplotypeCombinations = new ArrayList<>(haplotypeCombinations);
        this.wildTypeHaplotypeName = wildTypeHaplotypeName;
    }

    @Override
    public String toString()
    {
        return "HaplotypeAnalysis(" +
                "eventIdToCount=" + eventIdToCount +
                ", haplotypeCombinations=" + haplotypeCombinations +
                ')';
    }

    public List<String> getEventIds()
    {
        return new ArrayList<>(eventIdToCount.keySet());
    }

    public Optional<HaplotypeCombination> getBestHaplotypeCombination()
    {
        if (haplotypeCombinations.isEmpty())
            return Optional.empty();

        int minimumNowWildTypeCount = haplotypeCombinations.stream()
                .mapToInt(c -> c.getNonWildTypeCount(wildTypeHaplotypeName))
                .min().getAsInt();
        List<HaplotypeCombination> minimumCombinations = haplotypeCombinations.stream()
                .filter(c -> c.getNonWildTypeCount(wildTypeHaplotypeName) == minimumNowWildTypeCount)
                .collect(Collectors.toList());

        if (minimumCombinations.size() == 1)
            return Optional.of(minimumCombinations.get(0));
        else
            return Optional.empty();
    }

    public int getEventCount(String eventId)
    {
        return eventIdToCount.get(eventId);
    }
}
