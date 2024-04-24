package com.hartwig.hmftools.peach;

import static com.hartwig.hmftools.peach.PeachUtils.GERMLINE_TOTAL_COPY_NUMBER;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.peach.event.HaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;

import org.jetbrains.annotations.NotNull;

public class HaplotypeCaller
{
    public static @NotNull List<HaplotypeCombination> findPossibleHaplotypeCombinations(@NotNull Map<String, Integer> eventIdToCount,
            @NotNull List<NonDefaultHaplotype> nonDefaultHaplotypes, @NotNull DefaultHaplotype defaultHaplotype)
    {
        List<List<NonDefaultHaplotype>> nonDefaultCombinations =
                getPossibleNonDefaultHaplotypes(eventIdToCount, List.copyOf(nonDefaultHaplotypes));
        return nonDefaultCombinations.stream()
                .map(l -> constructHaplotypeCombination(l, defaultHaplotype))
                .collect(Collectors.toList());
    }

    @NotNull
    private static List<List<NonDefaultHaplotype>> getPossibleNonDefaultHaplotypes(@NotNull Map<String, Integer> eventIdToUnexplainedCount,
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
    private static HaplotypeCombination constructHaplotypeCombination(@NotNull List<NonDefaultHaplotype> nonDefaultHaplotypes,
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
}
