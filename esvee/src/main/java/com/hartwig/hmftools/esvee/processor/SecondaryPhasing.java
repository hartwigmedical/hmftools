package com.hartwig.hmftools.esvee.processor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;

public enum SecondaryPhasing
{
    ;

    public static <T extends SupportedAssembly> List<Set<T>> run(final List<Set<T>> primaryPhaseGroups)
    {
        // PERF: The O-notation on this sucks.
        final List<Set<T>> secondaryPhaseGroups = new ArrayList<>();
        for(Set<T> primaryPhaseGroup : primaryPhaseGroups)
        {
            if(primaryPhaseGroup.size() == 1)
            {
                secondaryPhaseGroups.add(primaryPhaseGroup);
                continue;
            }

            final Map<String, Set<T>> phasingByFragment = new HashMap<>();

            for(T assembly : primaryPhaseGroup)
                for(String fragment : assembly.getSupportFragments())
                    phasingByFragment.computeIfAbsent(fragment, __ -> new HashSet<>()).add(assembly);
            final Set<Set<T>> groups = phasingByFragment.values().stream()
                    .filter(g -> g.size() > 1)
                    .collect(Collectors.toSet());
            final Set<Set<T>> toRemove = new HashSet<>();
            for(Set<T> group : groups)
            {
                for(Set<T> other : groups)
                {
                    if(group == other)
                        continue;
                    if(other.containsAll(group))
                        toRemove.add(group);
                }
            }

            groups.stream()
                    .filter(g -> !toRemove.contains(g))
                    .forEach(secondaryPhaseGroups::add);
        }
        return secondaryPhaseGroups;
    }
}
