package com.hartwig.hmftools.esvee.processor;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.models.SupportedAssembly;

import org.jetbrains.annotations.Nullable;

public enum PrimaryPhasing
{
    ;

    public static <T extends SupportedAssembly> List<Set<T>> run(final List<T> assemblies)
    {
        final Map<String, Integer> primaryPhasingByFragment = new HashMap<>();
        final Map<Integer, Set<T>> primaryPhasing = new HashMap<>();
        int nextPhase = 0;

        for(final T assembly : assemblies)
        {
            int phasing = -1;
            final Set<Integer> additionalPhases = new HashSet<>();
            for(final String fragment : assembly.getSupportFragments())
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
                    for (final T otherAssembly : other)
                        for (final String otherFragment : otherAssembly.getSupportFragments())
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
            for(final String fragment : assembly.getSupportFragments())
                primaryPhasingByFragment.put(fragment, phasing);
        }

        return primaryPhasing.values().stream()
                .distinct()
                .collect(Collectors.toList());
    }
}
