package com.hartwig.hmftools.peach;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.NonWildTypeHaplotype;
import com.hartwig.hmftools.peach.haplotype.WildTypeHaplotype;
import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.stream.Collectors;

import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;
import static com.hartwig.hmftools.peach.PeachUtils.PCH_LOGGER;

public class HaplotypeCaller
{
    @NotNull
    private final HaplotypePanel haplotypePanel;

    public HaplotypeCaller(@NotNull HaplotypePanel haplotypePanel)
    {
        this.haplotypePanel = haplotypePanel;
    }

    public Map<String, List<HaplotypeCombination>> getGeneToPossibleHaplotypes(
            @NotNull Map<String, Integer> eventIdToCount
    )
    {
        Optional<String> nonPositiveCountEvent = eventIdToCount.entrySet().stream()
                .filter(e -> e.getValue() <= 0)
                .map(Map.Entry::getKey)
                .findAny();
        if (nonPositiveCountEvent.isPresent())
        {
            String nonPositiveCountEventName = nonPositiveCountEvent.get();
            String errorMsg = String.format(
                    "Events cannot have a non-positive count: %s -> %d",
                    nonPositiveCountEventName,
                    eventIdToCount.get(nonPositiveCountEventName)
            );
            throw new IllegalArgumentException(errorMsg);
        }
        return haplotypePanel.getGenes().stream()
                .collect(Collectors.toMap(g -> g, g -> getPossibleHaplotypes(eventIdToCount, g)));
    }

    private List<HaplotypeCombination> getPossibleHaplotypes(Map<String, Integer> eventIdToCount, String gene)
    {
        PCH_LOGGER.info("handling gene: {}", gene);
        Map<String, Integer> relevantEventIdToCount = eventIdToCount.entrySet().stream()
                .filter(e -> haplotypePanel.isRelevantFor(e.getKey(), gene))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
        PCH_LOGGER.info("events for gene '{}': {}", gene, relevantEventIdToCount);

        List<List<NonWildTypeHaplotype>> nonWildHaplotypeCombinations = getPossibleNonWildTypeHaplotypes(
                relevantEventIdToCount,
                List.copyOf(haplotypePanel.getNonWildTypeHaplotypes(gene))
        );
        return nonWildHaplotypeCombinations.stream()
                .map(l -> getCombination(l, haplotypePanel.getWildTypeHaplotype(gene)))
                .collect(Collectors.toList());
    }

    private List<List<NonWildTypeHaplotype>> getPossibleNonWildTypeHaplotypes(
            Map<String, Integer> eventIdToUnexplainedCount,
            List<NonWildTypeHaplotype> candidateHaplotypes)
    {
        // Use recursive descent to efficiently go through all possibilities
        if (eventIdToUnexplainedCount.values().stream().anyMatch(c -> c < 0))
        {
            return new ArrayList<>();
        }
        else if (eventIdToUnexplainedCount.values().stream().allMatch(c -> c == 0))
        {
            // return list containing only an empty list
            return new ArrayList<>(List.of(new ArrayList<>()));
        }

        List<List<NonWildTypeHaplotype>> results = new ArrayList<>();
        for (int i = 0; i < candidateHaplotypes.size(); i++)
        {
            NonWildTypeHaplotype haplotypeToTry = candidateHaplotypes.get(i);
            Map<String, Integer> eventIdToUnexplainedCountAfterTry = getEventIdToUnexplainedCountAfterTry(
                    eventIdToUnexplainedCount, haplotypeToTry
            );
            // To avoid encountering the exact same combination many times, limit the possible candidates for the recursive calls
            List<NonWildTypeHaplotype> candidateHaplotypesAfterTry = candidateHaplotypes.subList(i, candidateHaplotypes.size());
            List<List<NonWildTypeHaplotype>> subCombinations = getPossibleNonWildTypeHaplotypes(
                    eventIdToUnexplainedCountAfterTry,
                    candidateHaplotypesAfterTry
            );
            List<List<NonWildTypeHaplotype>> fullCombinations = subCombinations.stream()
                    .peek(l -> l.add(haplotypeToTry))
                    .collect(Collectors.toList());
            results.addAll(fullCombinations);
        }
        return results;
    }

    private Map<String, Integer> getEventIdToUnexplainedCountAfterTry(
            Map<String, Integer> eventIdToUnexplainedCount, NonWildTypeHaplotype haplotypeToTry
    )
    {
        HashMap<String, Integer> newEventIdToCount = Maps.newHashMap(eventIdToUnexplainedCount);
        for (HaplotypeEvent event : haplotypeToTry.events)
        {
            int newCount = newEventIdToCount.getOrDefault(event.id(), 0) - 1;
            newEventIdToCount.put(event.id(), newCount);
        }
        return newEventIdToCount;
    }

    private HaplotypeCombination getCombination(
            List<NonWildTypeHaplotype> nonWildTypeHaplotypes, WildTypeHaplotype wildTypeHaplotype
    )
    {
        Map<String, Integer> haplotypeNameToCount = nonWildTypeHaplotypes.stream()
                .map(h -> h.name)
                .collect(Collectors.groupingBy(h -> h, Collectors.summingInt(h -> 1)));
        if (nonWildTypeHaplotypes.size() < GERMLINE_TOTAL_COPY_NUMBER)
            haplotypeNameToCount.put(wildTypeHaplotype.name, GERMLINE_TOTAL_COPY_NUMBER - nonWildTypeHaplotypes.size());
        return new HaplotypeCombination(haplotypeNameToCount);
    }
}
