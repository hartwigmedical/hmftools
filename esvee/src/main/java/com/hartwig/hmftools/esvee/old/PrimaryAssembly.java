package com.hartwig.hmftools.esvee.old;

import static java.lang.String.format;

import com.hartwig.hmftools.esvee.common.Junction;

import org.jetbrains.annotations.Nullable;

public class PrimaryAssembly extends SupportedAssembly implements TrimmableAssembly<PrimaryAssembly>
{
    public final Junction OriginalJunction;
    public final String AnchorChromosome;
    public final int AnchorPosition;

    // if we started at this index and moved forward, we'd see the anchor
    public final int AnchorPositionInAssembly;

    public PrimaryAssembly(
            final String name, final String assembly, final Junction originalJunction,
            final String anchorChromosome, final int anchorPosition, final int anchorPositionInAssembly)
    {
        super(name, assembly);
        OriginalJunction = originalJunction;
        AnchorChromosome = anchorChromosome;
        AnchorPosition = anchorPosition;
        AnchorPositionInAssembly = anchorPositionInAssembly;
    }

    public PrimaryAssembly(
            final String name, final String assembly, final Junction originalJunction, final String anchorChromosome,
            final int anchorPosition, final int anchorPositionInAssembly, final PrimaryAssembly original)
    {
        this(name, assembly, originalJunction, anchorChromosome, anchorPosition, anchorPositionInAssembly);

        if(original.AnchorPosition == anchorPosition && original.AnchorChromosome.equals(anchorChromosome))
        {
            final int anchorDelta = anchorPositionInAssembly - original.AnchorPositionInAssembly;

            original.readSupport().forEach(x -> addEvidenceAt(x.Read, x.Index + anchorDelta));
        }
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

        for(ReadSupport support : readSupport())
        {
            final int newOffset = support.Index - removeLeft;
            if(newOffset >= newLength)
                continue;

            newAssembly.addEvidenceAt(support.Read, newOffset);
        }

        return newAssembly;
    }

    public String toString()
    {
        return format("%s junc(%s) support(%d)", Name, OriginalJunction, readSupportCount());
    }
}
