package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

import java.util.OptionalInt;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class SequenceUtils
{
    public static String buildSequence(final RefGenomeInterface refGenome, final SequenceDefinition definition)
    {
        String start = definition.startRegion() == null ? "" : getSequence(refGenome, definition.startRegion());
        if(definition.startOrientation() == Orientation.REVERSE)
        {
            start = reverseComplementBases(start);
        }
        String end = definition.endRegion() == null ? "" : getSequence(refGenome, definition.endRegion());
        if(definition.endOrientation() == Orientation.REVERSE)
        {
            end = reverseComplementBases(end);
        }
        return start + definition.insertSequence() + end;
    }

    private static String getSequence(final RefGenomeInterface refGenome, final ChrBaseRegion region)
    {
        String sequence = refGenome.getBaseString(region.chromosome(), region.start(), region.end());
        if(sequence == null || sequence.length() != region.baseLength())
        {
            throw new IllegalArgumentException("Attempt to create probe in unmapped region: " + region);
        }
        return sequence.toUpperCase();
    }

    public static boolean isDnaSequenceNormal(final String sequence)
    {
        return sequence.matches("^[acgtACGT]*$");
    }

    // Calculates the approximate size in bases of the insertion or deletion represented by the sequence.
    // Returns empty optional if it's a complex variant.
    public static OptionalInt sequenceIndelSize(final SequenceDefinition sequenceDefinition)
    {
        ChrBaseRegion start = sequenceDefinition.startRegion();
        ChrBaseRegion end = sequenceDefinition.endRegion();
        int insertLength = sequenceDefinition.insertSequence().length();

        if(start != null && end == null)
        {
            return OptionalInt.of(insertLength);
        }
        else if(start == null && end != null)
        {
            return OptionalInt.of(insertLength);
        }
        else if(start != null && end != null)
        {
            if(start.chromosome().equals(end.chromosome()))
            {
                // SNV, INDEL, or SV on same chromosome.
                if(start.start() > end.start())
                {
                    // Ensure start and end are ordered correctly to calculate the delete length.
                    start = sequenceDefinition.endRegion();
                    end = sequenceDefinition.startRegion();
                }
                // Clamp to >=0 because theoretically the regions could overlap in the case of an SV.
                int deleteLength = max(end.start() - start.end() - 1, 0);
                int difference = max(insertLength, deleteLength);
                return OptionalInt.of(difference);
            }
            else
            {
                // SV across different chromosomes.
                return OptionalInt.empty();
            }
        }
        else
        {
            // This shouldn't occur, but it would mean the probe is not based on the ref genome at all.
            return OptionalInt.empty();
        }
    }
}
