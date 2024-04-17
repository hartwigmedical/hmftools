package com.hartwig.hmftools.peach;

import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.function.Function;
import java.util.stream.Collector;
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

    @NotNull
    public Map<String, HaplotypeAnalysis> getGeneToHaplotypeAnalysis(@NotNull Map<String, Integer> eventIdToCount)
    {
        return haplotypePanel.getGenes().stream().collect(Collectors.toMap(g -> g, g -> getHaplotypeAnalysis(eventIdToCount, g)));
    }

    @NotNull
    private HaplotypeAnalysis getHaplotypeAnalysis(@NotNull Map<String, Integer> eventIdToCount, @NotNull String gene)
    {
        PCH_LOGGER.info("handle gene: {}", gene);
        Map<String, Integer> relevantEventIdToCount = eventIdToCount.entrySet()
                .stream()
                .filter(e -> haplotypePanel.isRelevantFor(e.getKey(), gene))
                .collect(toMapWithNull(Map.Entry::getKey, Map.Entry::getValue));

        List<List<NonDefaultHaplotype>> nonDefaultCombinations =
                getPossibleNonDefaultHaplotypes(relevantEventIdToCount, List.copyOf(haplotypePanel.getNonDefaultHaplotypes(gene)));
        List<HaplotypeCombination> possibleHaplotypeCombinations = nonDefaultCombinations.stream()
                .map(l -> getCombination(l, haplotypePanel.getDefaultHaplotype(gene)))
                .collect(Collectors.toList());
        return new HaplotypeAnalysis(relevantEventIdToCount, possibleHaplotypeCombinations, haplotypePanel.getDefaultHaplotype(gene)
                .getName(), haplotypePanel.getWildTypeHaplotypeName(gene));
    }

    @NotNull
    private List<List<NonDefaultHaplotype>> getPossibleNonDefaultHaplotypes(@NotNull Map<String, Integer> eventIdToUnexplainedCount,
            @NotNull List<NonDefaultHaplotype> candidateHaplotypes)
    {
        if(eventIdToUnexplainedCount.values().stream().anyMatch(Objects::isNull))
        {
            // cannot call haplotypes if count of a relevant event is unknown
            return new ArrayList<>();
        }

        // Use recursive descent to efficiently go through all possibilities
        if(eventIdToUnexplainedCount.values().stream().allMatch(c -> c == 0))
        {
            List<List<NonDefaultHaplotype>> possibleHaplotypeCombinations = new ArrayList<>();
            possibleHaplotypeCombinations.add(new ArrayList<>());
            return possibleHaplotypeCombinations;
        }

        assertNoNegativeEventCounts(eventIdToUnexplainedCount, candidateHaplotypes);

        List<List<NonDefaultHaplotype>> possibleHaplotypeCombinations = new ArrayList<>();
        for(int i = 0; i < candidateHaplotypes.size(); i++)
        {
            NonDefaultHaplotype haplotypeToTry = candidateHaplotypes.get(i);

            if(haplotypeIsPossible(eventIdToUnexplainedCount, haplotypeToTry))
            {
                Map<String, Integer> eventIdToUnexplainedCountAfterTry = eventIdToUnexplainedCount.keySet()
                        .stream()
                        .collect(Collectors.toMap(e -> e, e -> eventIdToUnexplainedCount.get(e) - haplotypeToTry.getMatchingCount(e)));

                // To avoid encountering the exact same combination many times, limit the possible candidates for the recursive calls
                List<NonDefaultHaplotype> candidateHaplotypesAfterTry = candidateHaplotypes.subList(i, candidateHaplotypes.size());

                List<List<NonDefaultHaplotype>> subCombinations =
                        getPossibleNonDefaultHaplotypes(eventIdToUnexplainedCountAfterTry, candidateHaplotypesAfterTry);
                List<List<NonDefaultHaplotype>> fullCombinations =
                        subCombinations.stream().peek(l -> l.add(haplotypeToTry)).collect(Collectors.toList());
                possibleHaplotypeCombinations.addAll(fullCombinations);
            }
        }
        return possibleHaplotypeCombinations;
    }

    private static void assertNoNegativeEventCounts(@NotNull Map<String, Integer> eventIdToUnexplainedCount,
            @NotNull List<NonDefaultHaplotype> candidateHaplotypes)
    {
        Optional<String> eventIdWithNegativeCount =
                eventIdToUnexplainedCount.entrySet().stream().filter(e -> e.getValue() < 0).map(Map.Entry::getKey).findFirst();
        if(eventIdWithNegativeCount.isPresent())
        {
            String candidateHaplotypesString =
                    candidateHaplotypes.stream().map(NonDefaultHaplotype::getName).collect(Collectors.joining(", "));
            throw new IllegalStateException(
                    "Negative count encountered for event '" + eventIdWithNegativeCount.get() + "' for candidate haplotypes ("
                            + candidateHaplotypesString + ")");
        }
    }

    private static boolean haplotypeIsPossible(@NotNull Map<String, Integer> eventIdToUnexplainedCount,
            @NotNull NonDefaultHaplotype haplotypeToTry)
    {
        return haplotypeToTry.events.stream()
                .map(HaplotypeEvent::id)
                .allMatch(e -> eventIdToUnexplainedCount.getOrDefault(e, 0) >= haplotypeToTry.getMatchingCount(e));
    }

    @NotNull
    private HaplotypeCombination getCombination(@NotNull List<NonDefaultHaplotype> nonDefaultHaplotypes,
            @NotNull DefaultHaplotype defaultHaplotype)
    {
        Map<String, Integer> haplotypeNameToCount =
                nonDefaultHaplotypes.stream().collect(Collectors.groupingBy(NonDefaultHaplotype::getName, Collectors.summingInt(h -> 1)));
        if(nonDefaultHaplotypes.size() < GERMLINE_TOTAL_COPY_NUMBER)
        {
            haplotypeNameToCount.put(defaultHaplotype.getName(), GERMLINE_TOTAL_COPY_NUMBER - nonDefaultHaplotypes.size());
        }
        return new HaplotypeCombination(haplotypeNameToCount);
    }

    @NotNull
    private static <T, K, U> Collector<T, ?, Map<K, U>> toMapWithNull(Function<? super T, ? extends K> keyMapper,
            Function<? super T, ? extends U> valueMapper)
    {
        // builtin toMap Collector cannot handle null values
        // source for this code: https://stackoverflow.com/questions/24630963/nullpointerexception-in-collectors-tomap-with-null-entry-values/32648397#32648397
        return Collectors.collectingAndThen(Collectors.toList(), list ->
        {
            Map<K, U> result = new HashMap<>();
            for(T item : list)
            {
                K key = keyMapper.apply(item);
                if(result.putIfAbsent(key, valueMapper.apply(item)) != null)
                {
                    throw new IllegalStateException(String.format("Duplicate key %s", key));
                }
            }
            return result;
        });
    }
}
