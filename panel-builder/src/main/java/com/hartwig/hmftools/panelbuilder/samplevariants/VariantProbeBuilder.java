package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.SequenceDefinition;

public class VariantProbeBuilder
{
    public static SequenceDefinition buildMutationProbe(final String chromosome, int position, final String ref, final String alt,
            int probeLength)
    {
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = probeLength / 2 - altLength / 2;

        ChrBaseRegion startRegion = regionEndingAt(chromosome, position - 1, startLength);
        int endBaseLength = probeLength - startRegion.baseLength() - altLength;

        int postPosition = position + refLength;
        ChrBaseRegion endRegion = regionStartingAt(chromosome, postPosition, endBaseLength);

        SequenceDefinition definition = SequenceDefinition.simpleVariant(startRegion, alt, endRegion);

        if(definition.baseLength() != probeLength)
        {
            throw new IllegalArgumentException(format("variant(%s:%d %s->%s) invalid sequenceLength(%d)",
                    chromosome, position, ref, alt, definition.baseLength()));
        }

        return definition;
    }

    public static SequenceDefinition buildSvProbe(final String chrStart, int positionStart, byte orientStart, final String chrEnd,
            int positionEnd, byte orientEnd, final String insertSequence, int probeLength)
    {
        int halfProbeLength = probeLength / 2;
        int insSeqLength = insertSequence.length();
        int halfInsSeqLength = insSeqLength / 2;
        int halfNonInsSeqLength = halfProbeLength - halfInsSeqLength;

        // +1 to -1 - do nothing
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        ChrBaseRegion startRegion;
        Orientation startOrient;
        ChrBaseRegion endRegion;
        Orientation endOrient;

        if(orientStart == ORIENT_FWD)
        {
            startRegion = regionEndingAt(chrStart, positionStart, halfNonInsSeqLength);
            startOrient = Orientation.FORWARD;

            int endBaseLength = probeLength - startRegion.baseLength() - insSeqLength;

            if(orientEnd == ORIENT_REV)
            {
                endRegion = regionStartingAt(chrEnd, positionEnd, endBaseLength);
                endOrient = Orientation.FORWARD;
            }
            else
            {
                endRegion = regionEndingAt(chrEnd, positionEnd, endBaseLength);
                endOrient = Orientation.REVERSE;
            }
        }
        else
        {
            if(orientEnd == ORIENT_FWD)
            {
                // swap ends and treat as +1/-1
                startRegion = regionEndingAt(chrEnd, positionEnd, halfNonInsSeqLength);
                startOrient = Orientation.FORWARD;

                int endBaseLength = probeLength - startRegion.baseLength() - insSeqLength;
                endRegion = regionStartingAt(chrStart, positionStart, endBaseLength);
            }
            else
            {
                // -1/-1 - start with the reversed bases from the end breakend
                startRegion = regionStartingAt(chrEnd, positionEnd, halfNonInsSeqLength);
                startOrient = Orientation.REVERSE;

                int endBaseLength = probeLength - startRegion.baseLength() - insSeqLength;
                endRegion = regionStartingAt(chrStart, positionStart, endBaseLength);
            }
            endOrient = Orientation.FORWARD;
        }

        SequenceDefinition definition =
                SequenceDefinition.structuralVariant(startRegion, startOrient, insertSequence, endRegion, endOrient);

        if(definition.baseLength() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return definition;
    }

    public static SequenceDefinition buildSglProbe(final String chromosome, int position, byte orientation, final String insertSequence,
            int probeLength)
    {
        int halfProbeLength = probeLength / 2;
        int insSeqLength = min(insertSequence.length(), halfProbeLength);
        int refBaseLength = probeLength - insSeqLength;

        ChrBaseRegion startRegion = null;
        ChrBaseRegion endRegion = null;
        String insert;

        if(orientation == ORIENT_FWD)
        {
            startRegion = regionEndingAt(chromosome, position, refBaseLength);
            insert = insertSequence.substring(0, insSeqLength);
        }
        else
        {
            insert = insertSequence.substring(insertSequence.length() - insSeqLength);
            endRegion = regionStartingAt(chromosome, position, refBaseLength);
        }

        SequenceDefinition definition = SequenceDefinition.simpleVariant(startRegion, insert, endRegion);

        if(definition.baseLength() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return definition;
    }
}
