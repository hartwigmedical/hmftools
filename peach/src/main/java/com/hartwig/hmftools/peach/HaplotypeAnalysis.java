package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class HaplotypeAnalysis
{
    @NotNull
    private final Map<String, Integer> eventIdToCount;
    @NotNull
    private final List<HaplotypeCombination> haplotypeCombinations;
    @NotNull
    private final String defaultHaplotypeName;
    @NotNull
    private final String wildTypeHaplotypeName;

    public HaplotypeAnalysis(
            @NotNull Map<String, Integer> eventIdToCount,
            @NotNull List<HaplotypeCombination> haplotypeCombinations,
            @NotNull String defaultHaplotypeName,
            @NotNull String wildTypeHaplotypeName
    )
    {
        this.eventIdToCount = new HashMap<>(eventIdToCount);
        this.haplotypeCombinations = new ArrayList<>(haplotypeCombinations);
        this.defaultHaplotypeName = defaultHaplotypeName;
        this.wildTypeHaplotypeName = wildTypeHaplotypeName;
    }

    @Override
    public String toString()
    {
        return "HaplotypeAnalysis(" +
                "eventIdToCount=" + eventIdToCount +
                ", haplotypeCombinations=" + haplotypeCombinations +
                ", defaultHaplotypeName=" + defaultHaplotypeName +
                ')';
    }

    public List<String> getEventIds()
    {
        return new ArrayList<>(eventIdToCount.keySet());
    }

    public List<HaplotypeCombination> getHaplotypeCombinations()
    {
        return new ArrayList<>(haplotypeCombinations);
    }

    @NotNull
    public String getDefaultHaplotypeName()
    {
        return defaultHaplotypeName;
    }

    public PeachQCStatus getAnalysisStatus()
    {
        if (haplotypeCombinations.isEmpty())
            return PeachQCStatus.FAIL_NO_COMBINATION_FOUND;

        List<HaplotypeCombination> minimumCombinations = getMinimumCombinations();
        if (minimumCombinations.size() != 1)
            return PeachQCStatus.FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND;
        else if (minimumCombinations.get(0).getHaplotypeCount() > 2)
            return PeachQCStatus.WARN_TOO_MANY_ALLELES_FOUND;
        else
            return PeachQCStatus.PASS;
    }

    public boolean hasBestHaplotypeCombination()
    {
        PeachQCStatus status = getAnalysisStatus();
        switch (status)
        {
            case PASS:
            case WARN_TOO_MANY_ALLELES_FOUND:
                return true;
            case FAIL_NO_COMBINATION_FOUND:
            case FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND:
                return false;
            default:
                throw new RuntimeException(String.format("Unrecognized QC status encountered: %s", status));
        }
    }

    public HaplotypeCombination getBestHaplotypeCombination()
    {
        if (!hasBestHaplotypeCombination())
            throw new RuntimeException(String.format("Gene does not have a best haplotype combination: status=%s", getAnalysisStatus()));

        return getMinimumCombinations().get(0);
    }

    @NotNull
    private List<HaplotypeCombination> getMinimumCombinations()
    {
        if (haplotypeCombinations.isEmpty())
            return Collections.emptyList();

        int minimumNonWildTypeCount = haplotypeCombinations.stream()
                .mapToInt(c -> c.getHaplotypeCountWithout(wildTypeHaplotypeName))
                .min().getAsInt();
        return haplotypeCombinations.stream()
                .filter(c -> c.getHaplotypeCountWithout(wildTypeHaplotypeName) == minimumNonWildTypeCount)
                .collect(Collectors.toList());
    }

    public int getEventCount(String eventId)
    {
        return eventIdToCount.get(eventId);
    }
}
