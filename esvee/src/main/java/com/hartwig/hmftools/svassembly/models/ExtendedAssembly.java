package com.hartwig.hmftools.svassembly.models;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.util.SequenceUtil;

public class ExtendedAssembly extends SupportedAssembly implements TrimmableAssembly<ExtendedAssembly>
{
    public final List<DiagramSet> Diagrams = new ArrayList<>();
    public final SupportedAssembly Source;

    public ExtendedAssembly(final String name, final String assembly, final SupportedAssembly source)
    {
        super(name, assembly);
        Source = source;
    }

    public void addDiagrams(@Nullable final DiagramSet diagrams)
    {
        if (diagrams != null)
            Diagrams.add(diagrams);
    }

    @Override
    public List<DiagramSet> getDiagrams()
    {
        final List<DiagramSet> diagrams = new ArrayList<>();
        diagrams.addAll(Source.getDiagrams());
        diagrams.addAll(Diagrams);
        return diagrams;
    }

    @Nullable
    @Override
    public ExtendedAssembly trim(final int removeLeft, final int removeRight)
    {
        final int newLength = getLength() - removeLeft - removeRight;
        if (newLength <= 0)
            return null;

        final String newBases = Assembly.substring(removeLeft, removeLeft + newLength);

        final ExtendedAssembly newAssembly = new ExtendedAssembly(Name, newBases, Source);
        newAssembly.Diagrams.addAll(Diagrams);
        for(final Map.Entry<Record, Integer> entry : getSupport())
        {
            final int newOffset = entry.getValue() - removeLeft;
            if (newOffset >= newLength)
                continue;

            newAssembly.addEvidenceAt(entry.getKey(), newOffset);
        }

        return newAssembly;
    }

    public ExtendedAssembly flipStrand()
    {
        final String assembly = SequenceUtil.reverseComplement(Assembly);
        final ExtendedAssembly flipped = new ExtendedAssembly(Name, assembly, Source);

        for (final Map.Entry<Record, Integer> support : getSupport())
            flipped.addEvidenceAt(support.getKey().flipStrand(),
                    getLength() - support.getValue() - support.getKey().getLength());
        flipped.recalculateBaseQuality();
        return flipped;
    }
}
