package com.hartwig.hmftools.esvee.processor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.sequence.ExtendedAssembly;

public final class SecondaryPhasing
{
    public static List<Set<ExtendedAssembly>> run(final List<Set<ExtendedAssembly>> primaryPhaseGroups)
    {
        // PERF: The O-notation on this sucks.
        final List<Set<ExtendedAssembly>> secondaryPhaseGroups = new ArrayList<>();
        for(Set<ExtendedAssembly> primaryPhaseGroup : primaryPhaseGroups)
        {
            if(primaryPhaseGroup.size() == 1)
            {
                secondaryPhaseGroups.add(primaryPhaseGroup);
                continue;
            }

            final Map<String, Set<ExtendedAssembly>> phasingByFragment = new HashMap<>();

            for(ExtendedAssembly assembly : primaryPhaseGroup)
            {
                for(String fragment : assembly.getSupportReadNames())
                {
                    phasingByFragment.computeIfAbsent(fragment, __ -> new HashSet<>()).add(assembly);
                }
            }

            Set<Set<ExtendedAssembly>> groups = phasingByFragment.values().stream()
                    .filter(g -> g.size() > 1)
                    .collect(Collectors.toSet());

            Set<Set<ExtendedAssembly>> toRemove = new HashSet<>();

            for(Set<ExtendedAssembly> group : groups)
            {
                for(Set<ExtendedAssembly> other : groups)
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
