package com.hartwig.hmftools.esvee.old;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.jetbrains.annotations.Nullable;

public final class PrimaryPhasing
{
    public static <T extends SupportedAssembly> List<Set<T>> run(final List<T> assemblies)
    {
        final Map<String, Integer> primaryPhasingByFragment = new HashMap<>();
        final Map<Integer, Set<T>> primaryPhasing = new HashMap<>();
        int nextPhase = 0;

        for(T assembly : assemblies)
        {
            int phasing = -1;
            final Set<Integer> additionalPhases = new HashSet<>();
            for(String fragment : assembly.getSupportReadNames())
            {
                @Nullable
                final Integer fragmentPhase = primaryPhasingByFragment.get(fragment);
                if(fragmentPhase == null || fragmentPhase == phasing)
                    continue;

                if(phasing == -1)
                    phasing = fragmentPhase;
                else if(additionalPhases.add(fragmentPhase))
                {
                    final Set<T> current = primaryPhasing.get(phasing);
                    final Set<T> other = primaryPhasing.get(fragmentPhase);

                    if(current == other)
                        continue;
                    current.addAll(other);
                    for(T otherAssembly : other)
                        for(String otherFragment : otherAssembly.getSupportReadNames())
                            primaryPhasingByFragment.put(otherFragment, phasing);
                    primaryPhasing.remove(fragmentPhase);
                }
            }

            if(phasing == -1)
            {
                phasing = nextPhase++;
                primaryPhasing.put(phasing, Collections.newSetFromMap(new IdentityHashMap<>()));
            }
            primaryPhasing.get(phasing).add(assembly);
            for(String fragment : assembly.getSupportReadNames())
                primaryPhasingByFragment.put(fragment, phasing);
        }

        return primaryPhasing.values().stream()
                .distinct()
                .collect(Collectors.toList());
    }
}
