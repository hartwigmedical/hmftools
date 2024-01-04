package com.hartwig.hmftools.esvee.assembly;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.ToIntFunction;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.sequence.AlignedAssembly;
import com.hartwig.hmftools.esvee.sequence.Alignment;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Builds the anchor cigar for a variant
public class AnchorCigarFactory
{
    public AnchorCigarFactory()
    {
    }

    public Pair<Pair<Cigar, Integer>, Pair<Cigar, Integer>> anchorCigar(final AlignedAssembly assembly,
            @Nullable final Alignment left, @Nullable final Alignment right)
    {
        return anchorCigar(assembly.getAlignmentBlocks(), left, right);
    }

    @VisibleForTesting
    Pair<Pair<Cigar, Integer>, Pair<Cigar, Integer>> anchorCigar(final List<Alignment> assemblyAlignment,
            @Nullable final Alignment left, @Nullable final Alignment right)
    {
        return Pair.of(anchorCigarLeft(assemblyAlignment, left), anchorCigarRight(assemblyAlignment, right));
    }

    @Nullable
    private Pair<Cigar, Integer> anchorCigarLeft(final List<Alignment> assemblyAlignment, @Nullable final Alignment alignment)
    {
        return computeAnchorCigar(assemblyAlignment, alignment,
                -1,
                a -> a.Inverted ? a.ReferenceStartPosition : a.ReferenceStartPosition + a.Length - 1,
                a -> a.Inverted ? a.ReferenceStartPosition + a.Length - 1 : a.ReferenceStartPosition);
    }

    @Nullable
    private Pair<Cigar, Integer> anchorCigarRight(final List<Alignment> assemblyAlignment, @Nullable final Alignment alignment)
    {
        return computeAnchorCigar(assemblyAlignment, alignment,
                1, a -> a.ReferenceStartPosition, a -> a.ReferenceStartPosition + a.Length - 1);
    }

    @Nullable
    private Pair<Cigar, Integer> computeAnchorCigar(final List<Alignment> assemblyAlignment, @Nullable final Alignment alignment,
            final int step, final ToIntFunction<Alignment> startComputer, final ToIntFunction<Alignment> endComputer)
    {
        if(alignment == null)
            return null;

        final boolean inverted = alignment.Inverted;

        final List<CigarElement> elements = new ArrayList<>();
        elements.add(new CigarElement(alignment.Length, CigarOperator.M));
        // Where are we at the moment
        int referenceIndex = endComputer.applyAsInt(alignment);

        final int startIndex = assemblyAlignment.indexOf(alignment);
        int index = startIndex + step;
        for(; index < assemblyAlignment.size() && index >= 0; index += step)
        {
            final Alignment current = assemblyAlignment.get(index);

            if(current.isUnmapped())
            {
                if(current.Length >= SvConstants.VARIANT_MIN_LENGTH)
                    break;
                elements.add(new CigarElement(current.Length, CigarOperator.I));
                continue;
            }
            if(current.Inverted != inverted)
                break;

            final int skip = (startComputer.applyAsInt(current) - referenceIndex) * step - 1;
            if(skip < 0)
                break; // Duplication
            else if(skip >= SvConstants.VARIANT_MIN_LENGTH)
                break; // Structural variant
            else if(skip > 0)
                elements.add(new CigarElement(skip, CigarOperator.D));

            elements.add(new CigarElement(current.Length, CigarOperator.M));
            referenceIndex = endComputer.applyAsInt(current);
        }

        // Now compute leftover
        int leftover = 0;
        while(elements.get(elements.size() - 1).getOperator() == CigarOperator.I)
        {
            leftover += elements.get(elements.size() - 1).getLength();
            elements.remove(elements.size() - 1);
        }
        for(; index < assemblyAlignment.size() && index >= 0; index += step)
            leftover += assemblyAlignment.get(index).Length;
        if(leftover > 0)
            elements.add(new CigarElement(leftover, CigarOperator.S));

        boolean reverseCigar = step == -1;
        if(alignment.Inverted)
            reverseCigar = !reverseCigar;
        if(reverseCigar)
            Collections.reverse(elements);

        final Cigar cigar = new Cigar(elements);
        return Pair.of(cigar, Cigar.getReadLength(elements));
    }
}
