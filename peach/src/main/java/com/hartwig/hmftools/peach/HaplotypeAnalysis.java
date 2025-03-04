package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.jetbrains.annotations.Nullable;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class HaplotypeAnalysis
{
    private final Map<String, Integer> eventIdToCount;
    private final List<HaplotypeCombination> haplotypeCombinations;
    private final String defaultHaplotypeName;
    private final String wildTypeHaplotypeName;
    private final PeachQCStatus qcStatus;
    private final HaplotypeCombination bestHaplotypeCombination;

    public HaplotypeAnalysis(final Map<String, Integer> eventIdToCount, final List<HaplotypeCombination> haplotypeCombinations,
            final String defaultHaplotypeName, final String wildTypeHaplotypeName, final PeachQCStatus qcStatus,
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
    public String toString()
    {
        return "HaplotypeAnalysis(" + "eventIdToCount=" + eventIdToCount + ", haplotypeCombinations=" + haplotypeCombinations
                + ", defaultHaplotypeName=" + defaultHaplotypeName + ')';
    }

    public Set<String> getEventIds()
    {
        return new HashSet<>(eventIdToCount.keySet());
    }

    public List<HaplotypeCombination> getHaplotypeCombinations()
    {
        return new ArrayList<>(haplotypeCombinations);
    }

    public String getWildTypeHaplotypeName()
    {
        return wildTypeHaplotypeName;
    }

    public String getDefaultHaplotypeName()
    {
        return defaultHaplotypeName;
    }

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
    public Integer getEventCount(final String eventId)
    {
        return eventIdToCount.get(eventId);
    }
}
