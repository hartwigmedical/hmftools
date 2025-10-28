package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;

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
        String insert = definition.insertSequence() == null ? "" : definition.insertSequence();
        String end = definition.endRegion() == null ? "" : getSequence(refGenome, definition.endRegion());
        if(definition.endOrientation() == Orientation.REVERSE)
        {
            end = reverseComplementBases(end);
        }
        return start + insert + end;
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
}
