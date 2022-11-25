package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class HaplotypeCaller
{
    @NotNull
    private final HaplotypePanel haplotypePanel;

    public HaplotypeCaller(
            @NotNull HaplotypePanel haplotypePanel
    )
    {
        this.haplotypePanel = haplotypePanel;
    }

    public void callPossibleHaplotypes(Map<HaplotypeEvent, Integer> eventToCount)
    {
        //haplotypePanel.getGenes().stream().map(g -> callPossibleHaplotypes(eventToCount, g));

        Map<HaplotypeEvent, Set<String>> eventToRelevantGenes = getEventToRelevantGenes(eventToCount.keySet());
    }

    private void callPossibleHaplotypes(Map<HaplotypeEvent, Integer> eventToCount, String gene)
    {
        PCH_LOGGER.info("handling gene: {}", gene);
        Map<HaplotypeEvent, Integer> relevantEventToCount = eventToCount.entrySet().stream()
                .filter(e -> haplotypePanel.isRelevantFor(e.getKey(), gene))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        PCH_LOGGER.info("events for gene '{}': {}", gene, relevantEventToCount);

    }

    @NotNull
    private Map<HaplotypeEvent, Set<String>> getEventToRelevantGenes(Set<HaplotypeEvent> events)
    {
        return events.stream().collect(Collectors.toMap(e -> e, haplotypePanel::getRelevantGenes));
    }


}
