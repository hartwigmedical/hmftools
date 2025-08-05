package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class VariantProbeBuilder
{
    public static VariantProbeData buildMutationProbe(final String chromosome, final int position, final String ref, final String alt,
            final RefGenomeInterface refGenome)
    {
        int probeLength = PROBE_LENGTH;
        int altLength = alt.length();
        int refLength = ref.length();
        int startLength = probeLength / 2 - altLength / 2;
        int startPos = position - startLength;

        ChrBaseRegion startRegion = new ChrBaseRegion(chromosome, startPos, position - 1);
        String startBases = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
        int endBaseLength = probeLength - startBases.length() - altLength;

        int postPosition = position + refLength;
        ChrBaseRegion endRegion = new ChrBaseRegion(chromosome, postPosition, postPosition + endBaseLength - 1);
        String endBases = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());

        String sequence = startBases + alt + endBases;

        if(sequence.length() != probeLength)
        {
            throw new IllegalArgumentException(format("variant(%s:%d %s->%s) invalid sequenceLength(%d): %s",
                    chromosome, position, ref, alt, sequence.length(), sequence));
        }

        return new VariantProbeData(sequence, startRegion, alt, endRegion);
    }

    public static VariantProbeData buildSvProbe(final String chrStart, final int positionStart, final byte orientStart,
            final String chrEnd, final int positionEnd, final byte orientEnd, final String insertSequence,
            final RefGenomeInterface refGenome)
    {
        int probeLength = PROBE_LENGTH;
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
            int probeStart = positionStart - halfNonInsSeqLength + 1;
            startRegion = new ChrBaseRegion(chrStart, probeStart, positionStart);
            basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());

            int endBaseLength = probeLength - basesStart.length() - insSeqLength;

            if(orientEnd == ORIENT_REV)
            {
                endRegion = new ChrBaseRegion(chrEnd, positionEnd, positionEnd + endBaseLength - 1);
                basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
            }
            else
            {
                endRegion = new ChrBaseRegion(chrEnd, positionEnd - endBaseLength + 1, positionEnd);
                basesEnd = refGenome.getBaseString(endRegion.chromosome(), endRegion.start(), endRegion.end());
                basesEnd = Nucleotides.reverseComplementBases(basesEnd);
            }
        }
        else
        {
            if(orientEnd == ORIENT_FWD)
            {
                // swap ends and treat as +1/-1
                int probeStart = positionEnd - halfNonInsSeqLength + 1;
                startRegion = new ChrBaseRegion(chrEnd, probeStart, positionEnd);
                basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                endRegion = new ChrBaseRegion(chrStart, positionStart, positionStart + endBaseLength - 1);
            }
            else
            {
                // -1/-1 - start with the reversed bases from the end breakend
                startRegion = new ChrBaseRegion(chrEnd, positionEnd, positionEnd + halfNonInsSeqLength - 1);
                basesStart = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
                basesStart = Nucleotides.reverseComplementBases(basesStart);

                int endBaseLength = probeLength - basesStart.length() - insSeqLength;
                endRegion = new ChrBaseRegion(chrStart, positionStart, positionStart + endBaseLength - 1);
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

    public static VariantProbeData buildSglProbe(final String chromosome, final int position, final byte orientation,
            final String insertSequence, final RefGenomeInterface refGenome)
    {
        int probeLength = PROBE_LENGTH;
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
            int probeStart = position - refBaseLength + 1;
            startRegion = new ChrBaseRegion(chromosome, probeStart, position);
            String refBases = refGenome.getBaseString(startRegion.chromosome(), startRegion.start(), startRegion.end());
            insert = insertSequence.substring(0, insSeqLength);
            sequence = refBases + insert;
        }
        else
        {
            insert = insertSequence.substring(insertSequence.length() - insSeqLength);
            endRegion = new ChrBaseRegion(chromosome, position, position + refBaseLength - 1);
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
