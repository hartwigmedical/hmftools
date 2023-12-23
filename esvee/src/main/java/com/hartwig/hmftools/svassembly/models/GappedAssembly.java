package com.hartwig.hmftools.svassembly.models;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.svassembly.assembly.SupportChecker;

import org.jetbrains.annotations.Nullable;

public class GappedAssembly extends SupportedAssembly
{
    public final List<ExtendedAssembly> Sources;

    public GappedAssembly(final String name, final List<ExtendedAssembly> sources)
    {
        super(name, createMergedAssembly(sources));
        Sources = sources;
    }

    @Override
    public boolean tryAddSupport(final SupportChecker checker, final Record record)
    {
        int currentOffset = 0;
        for (final ExtendedAssembly assembly : Sources)
        {
            @Nullable
            final Integer index = checker.WeakSupport.supportIndex(assembly, record);

            if (index != null)
            {
                addEvidenceAt(record, currentOffset + index);
                return true;
            }

            currentOffset += assembly.Assembly.length() + 1;
        }

        return false;
    }

    private static String createMergedAssembly(final Collection<ExtendedAssembly> sources)
    {
        return sources.stream()
                .map(s -> s.Assembly)
                .collect(Collectors.joining("X"));
    }

    @Override
    public List<DiagramSet> getDiagrams()
    {
        return Sources.stream()
                .flatMap(assembly -> assembly.getDiagrams().stream())
                .collect(Collectors.toList());
    }

    public GappedAssembly flipStrand()
    {
        final List<ExtendedAssembly> newSources = Sources.stream()
                .map(ExtendedAssembly::flipStrand)
                .collect(Collectors.toList());
        Collections.reverse(newSources);

        final GappedAssembly flipped = new GappedAssembly(Name, newSources);
        for (final Map.Entry<Record, Integer> support : getSupport())
            flipped.addEvidenceAt(support.getKey().flipStrand(),
                    getLength() - support.getValue() - support.getKey().getLength());
        flipped.recalculateBaseQuality();
        return flipped;
    }
}
