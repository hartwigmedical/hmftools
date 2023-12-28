package com.hartwig.hmftools.esvee.models;

import static com.hartwig.hmftools.esvee.read.ReadUtils.flipRead;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.html.DiagramSet;
import com.hartwig.hmftools.esvee.read.Read;

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
        if(diagrams != null)
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
        if(newLength <= 0)
            return null;

        final String newBases = Assembly.substring(removeLeft, removeLeft + newLength);

        final ExtendedAssembly newAssembly = new ExtendedAssembly(Name, newBases, Source);
        newAssembly.Diagrams.addAll(Diagrams);
        for(Map.Entry<Read, Integer> entry : getSupport())
        {
            final int newOffset = entry.getValue() - removeLeft;
            if(newOffset >= newLength)
                continue;

            newAssembly.addEvidenceAt(entry.getKey(), newOffset);
        }

        return newAssembly;
    }

    public ExtendedAssembly flipStrand()
    {
        final String assembly = SequenceUtil.reverseComplement(Assembly);
        final ExtendedAssembly flipped = new ExtendedAssembly(Name, assembly, Source);

        for(Map.Entry<Read, Integer> support : getSupport())
        {
            int initialReadLength = support.getKey().getLength();
            Read flippedRead = flipRead(support.getKey());

            flipped.addEvidenceAt(flippedRead,getLength() - support.getValue() - initialReadLength);
        }
        flipped.recalculateBaseQuality();
        return flipped;
    }
}
