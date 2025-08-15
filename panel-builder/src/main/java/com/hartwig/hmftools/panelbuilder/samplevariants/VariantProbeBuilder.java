package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionEndingAt;
import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionStartingAt;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class VariantProbeBuilder
{
    public static VariantProbeData buildMutationProbe(final String chromosome, int position, final String ref, final String alt,
            int probeLength, final RefGenomeInterface refGenome)
    {
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = probeLength / 2 - altLength / 2;

        ChrBaseRegion startRegion = regionEndingAt(chromosome, position - 1, startLength);
        String startBases = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
        int endBaseLength = probeLength - startBases.length() - altLength;

        int postPosition = position + refLength;
        ChrBaseRegion endRegion = regionStartingAt(chromosome, postPosition, endBaseLength);
        String endBases = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());

        String sequence = startBases + alt + endBases;

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException(format("variant(%s:%d %s->%s) invalid sequenceLength(%d): %s",
                    chromosome, position, ref, alt, sequence.length(), sequence));
        }

        return new VariantProbeData(sequence, startRegion, alt, endRegion);
    }

    public static VariantProbeData buildSvProbe(final String chrStart, int positionStart, byte orientStart, final String chrEnd,
            int positionEnd, byte orientEnd, final String insertSequence, int probeLength, final RefGenomeInterface refGenome)
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
        String basesStart;
        ChrBaseRegion endRegion;
        String basesEnd;

        if(orientStart == ORIENT_FWD)
        {
            startRegion = regionEndingAt(chrStart, positionStart, halfNonInsSeqLength);
            basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());

            int endBaseLength = probeLength - basesStart.length() - insSeqLength;

            if(orientEnd == ORIENT_REV)
            {
                endRegion = regionStartingAt(chrEnd, positionEnd, endBaseLength);
                basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
            }
            else
            {
                endRegion = regionEndingAt(chrEnd, positionEnd, endBaseLength);
                basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
                basesEnd = reverseComplementBases(basesEnd);
            }
        }
        else
        {
            if(orientEnd == ORIENT_FWD)
            {
                // swap ends and treat as +1/-1
                startRegion = regionEndingAt(chrEnd, positionEnd, halfNonInsSeqLength);
                basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                endRegion = regionStartingAt(chrStart, positionStart, endBaseLength);
            }
            else
            {
                // -1/-1 - start with the reversed bases from the end breakend
                startRegion = regionStartingAt(chrEnd, positionEnd, halfNonInsSeqLength);
                basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
                basesStart = reverseComplementBases(basesStart);

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                endRegion = regionStartingAt(chrStart, positionStart, endBaseLength);
            }
            basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
        }

        String sequence = basesStart + insertSequence + basesEnd;

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return new VariantProbeData(sequence, startRegion, insertSequence, endRegion);
    }

    public static VariantProbeData buildSglProbe(final String chromosome, int position, byte orientation, final String insertSequence,
            int probeLength, final RefGenomeInterface refGenome)
    {
        int halfProbeLength = probeLength / 2;
        int insSeqLength = min(insertSequence.length(), halfProbeLength);
        int refBaseLength = probeLength - insSeqLength;

        // +1 take as-is
        // -1 to +1 - a DUP

        // +1 to +1 - start normal, add insert and then reverse compliment of the other side
        // -1 to -1 - alt insert sequence then ref

        String sequence;

        ChrBaseRegion startRegion = null;
        ChrBaseRegion endRegion = null;
        String insert;

        if(orientation == ORIENT_FWD)
        {
            startRegion = regionEndingAt(chromosome, position, refBaseLength);
            String refBases = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
            insert = insertSequence.substring(0, insSeqLength);
            sequence = refBases + insert;
        }
        else
        {
            insert = insertSequence.substring(insertSequence.length() - insSeqLength);
            endRegion = regionStartingAt(chromosome, position, refBaseLength);
            String refBases = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
            sequence = insert + refBases;
        }

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException("Invalid probe");
        }

        return new VariantProbeData(sequence, startRegion, insert, endRegion);
    }
}
