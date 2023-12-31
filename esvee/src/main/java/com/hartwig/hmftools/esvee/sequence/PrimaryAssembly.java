package com.hartwig.hmftools.esvee.sequence;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.html.DiagramSet;
import com.hartwig.hmftools.esvee.read.Read;

import org.jetbrains.annotations.Nullable;

public class PrimaryAssembly extends SupportedAssembly implements TrimmableAssembly<PrimaryAssembly>
{
    public final Junction OriginalJunction;
    public final String AnchorChromosome;
    public final int AnchorPosition;

    // if we started at this index and moved forward, we'd see the anchor
    public final int AnchorPositionInAssembly;

    public final List<DiagramSet> Diagrams;

    public PrimaryAssembly(
            final String name, final String assembly, final Junction originalJunction,
            final String anchorChromosome, final int anchorPosition, final int anchorPositionInAssembly)
    {
        super(name, assembly);
        OriginalJunction = originalJunction;
        AnchorChromosome = anchorChromosome;
        AnchorPosition = anchorPosition;
        AnchorPositionInAssembly = anchorPositionInAssembly;
        Diagrams = new ArrayList<>();
    }

    public PrimaryAssembly(
            final String name, final String assembly, final Junction originalJunction, final String anchorChromosome,
            final int anchorPosition, final int anchorPositionInAssembly, final PrimaryAssembly original)
    {
        this(name, assembly, originalJunction, anchorChromosome, anchorPosition, anchorPositionInAssembly);

        Diagrams.addAll(original.Diagrams);

        if(original.AnchorPosition == anchorPosition && original.AnchorChromosome.equals(anchorChromosome))
        {
            final int anchorDelta = anchorPositionInAssembly - original.AnchorPositionInAssembly;

            original.readSupport().forEach(x -> addEvidenceAt(x.Read, x.Index + anchorDelta));
        }
    }

    public void addDiagrams(@Nullable final DiagramSet diagrams)
    {
        if(diagrams != null)
            Diagrams.add(diagrams);
    }

    @Override
    public List<DiagramSet> getDiagrams()
    {
        return Diagrams;
    }

    @Nullable
    @Override
    public PrimaryAssembly trim(final int removeLeft, final int removeRight)
    {
        final int newLength = getLength() - removeLeft - removeRight;
        if(newLength <= 0)
            return null;

        final String newBases = Assembly.substring(removeLeft, removeLeft + newLength);

        int newAnchorPositionInAssembly = AnchorPositionInAssembly + removeLeft;
        int newAnchorPosition = AnchorPosition;
        if(newAnchorPositionInAssembly >= newBases.length())
        {
            // Move anchor left
            final int leftMove = newBases.length() - newAnchorPositionInAssembly + 1;
            newAnchorPositionInAssembly -= leftMove;
            newAnchorPosition -= leftMove;
        }

        final PrimaryAssembly newAssembly = new PrimaryAssembly(
                Name, newBases, OriginalJunction, AnchorChromosome, newAnchorPosition, newAnchorPositionInAssembly);

        newAssembly.Diagrams.addAll(Diagrams);

        for(ReadSupport support : readSupport())
        {
            final int newOffset = support.Index - removeLeft;
            if(newOffset >= newLength)
                continue;

            newAssembly.addEvidenceAt(support.Read, newOffset);
        }

        return newAssembly;
    }
}
