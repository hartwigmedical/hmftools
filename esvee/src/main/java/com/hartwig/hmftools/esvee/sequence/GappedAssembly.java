package com.hartwig.hmftools.esvee.sequence;

import static com.hartwig.hmftools.esvee.read.ReadUtils.flipRead;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.assembly.SupportChecker;
import com.hartwig.hmftools.esvee.html.DiagramSet;
import com.hartwig.hmftools.esvee.read.Read;

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
    public boolean tryAddSupport(final SupportChecker checker, final Read read)
    {
        int currentOffset = 0;
        for(ExtendedAssembly assembly : Sources)
        {
            @Nullable
            final Integer index = checker.WeakSupport.supportIndex(assembly, read);

            if(index != null)
            {
                addEvidenceAt(read, currentOffset + index);
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

        for(ReadSupport support : readSupport())
        {
            int initialReadLength = support.Read.getLength();
            Read flippedRead = flipRead(support.Read);
            flipped.addEvidenceAt(flippedRead,getLength() - support.Index - initialReadLength);
        }

        flipped.recalculateBaseQuality();
        return flipped;
    }
}
