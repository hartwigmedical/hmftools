package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.event.HaplotypeEventFactory;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class HaplotypeCaller
{
    @NotNull
    private final HaplotypePanel haplotypePanel;

    public HaplotypeCaller(@NotNull HaplotypePanel haplotypePanel)
    {
        this.haplotypePanel = haplotypePanel;
    }

    public void callPossibleHaplotypes(@NotNull Map<String, Integer> eventIdToCount)
    {
        //haplotypePanel.getGenes().stream().map(g -> callPossibleHaplotypes(eventIdToCount, g));
    }

    private void callPossibleHaplotypes(Map<String, Integer> eventIdToCount, String gene)
    {
        PCH_LOGGER.info("handling gene: {}", gene);
        Map<String, Integer> relevantEventIdToCount = eventIdToCount.entrySet().stream()
                .filter(e -> haplotypePanel.isRelevantFor(e.getKey(), gene))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        PCH_LOGGER.info("events for gene '{}': {}", gene, relevantEventIdToCount);
        
    }
}
