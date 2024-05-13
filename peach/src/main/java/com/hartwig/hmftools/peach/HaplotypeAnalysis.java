package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    @NotNull
    private final PeachQCStatus qcStatus;
    @Nullable
    private final HaplotypeCombination bestHaplotypeCombination;

    public HaplotypeAnalysis(@NotNull Map<String, Integer> eventIdToCount, @NotNull List<HaplotypeCombination> haplotypeCombinations,
            @NotNull String defaultHaplotypeName, @NotNull String wildTypeHaplotypeName, @NotNull PeachQCStatus qcStatus,
            @Nullable HaplotypeCombination bestHaplotypeCombination)
    {
        this.eventIdToCount = new HashMap<>(eventIdToCount);
        this.haplotypeCombinations = new ArrayList<>(haplotypeCombinations);
        this.defaultHaplotypeName = defaultHaplotypeName;
        this.wildTypeHaplotypeName = wildTypeHaplotypeName;
        this.qcStatus = qcStatus;
        this.bestHaplotypeCombination = bestHaplotypeCombination;
    }

    @Override
    @NotNull
    public String toString()
    {
        return "HaplotypeAnalysis(" + "eventIdToCount=" + eventIdToCount + ", haplotypeCombinations=" + haplotypeCombinations
                + ", defaultHaplotypeName=" + defaultHaplotypeName + ')';
    }

    @NotNull
    public Set<String> getEventIds()
    {
        return new HashSet<>(eventIdToCount.keySet());
    }

    @NotNull
    public List<HaplotypeCombination> getHaplotypeCombinations()
    {
        return new ArrayList<>(haplotypeCombinations);
    }

    @NotNull
    public String getWildTypeHaplotypeName()
    {
        return wildTypeHaplotypeName;
    }

    @NotNull
    public String getDefaultHaplotypeName()
    {
        return defaultHaplotypeName;
    }

    @NotNull
    public PeachQCStatus getQcStatus()
    {
        return qcStatus;
    }

    @Nullable
    public HaplotypeCombination getBestHaplotypeCombination()
    {
        return bestHaplotypeCombination;
    }

    @Nullable
    public Integer getEventCount(@NotNull String eventId)
    {
        return eventIdToCount.get(eventId);
    }
}
