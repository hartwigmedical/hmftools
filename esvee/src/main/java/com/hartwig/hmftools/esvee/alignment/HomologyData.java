package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.alignment.BreakendBuilder.segmentOrientation;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class HomologyData
{
    public final String Homology;
    public final int ExactStart;
    public final int ExactEnd;
    public final int InexactStart;
    public final int InexactEnd;

    public HomologyData(final String homology, final int exactStart, final int exactEnd, final int inexactStart, final int inexactEnd)
    {
        Homology = homology;
        ExactStart = exactStart;
        ExactEnd = exactEnd;
        InexactStart = inexactStart;
        InexactEnd = inexactEnd;
    }

    public String toString() { return format("%s exact(%d,%d) inexact(%d,%d)", Homology, ExactStart, ExactEnd, InexactStart, ExactEnd); }

    public static HomologyData determineHomology(
            final AlignData alignStart, final AlignData alignEnd, final RefGenomeInterface refGenome)
    {
        int overlap = alignStart.sequenceEnd() - alignEnd.sequenceStart() + 1;

        if(overlap <= 0)
            return null;

        byte orientationStart = segmentOrientation(alignStart, true);
        String basesStart = getOverlapBases(alignStart, orientationStart, overlap, refGenome);

        byte orientationEnd = segmentOrientation(alignEnd, false);
        String basesEnd = getOverlapBases(alignEnd, orientationEnd, overlap, refGenome);

        if(orientationStart == orientationEnd)
        {
            basesEnd = Nucleotides.reverseComplementBases(basesEnd);
        }

        return determineHomology(basesStart, basesEnd, overlap);
    }

    public static HomologyData determineHomology(final String basesStart, final String basesEnd, final int overlap)
    {
        int exactMatch = 0;
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < overlap; ++i)
        {
            if(basesStart.charAt(i) == basesEnd.charAt(i))
            {
                sb.append(basesStart.charAt(i));
                ++exactMatch;
            }
            else
            {
                break;
            }
        }

        if(exactMatch == 0)
            return new HomologyData("", 0, 0, 0, overlap);

        int exactMid = exactMatch / 2;
        int exactStart = exactMid;
        int exactEnd = exactMatch - exactStart;

        int inexactStart = 0;
        int inexactEnd = overlap - exactMatch;

        return new HomologyData(sb.toString(), -exactStart, exactEnd, -inexactStart, inexactEnd);
    }

    private static String getOverlapBases(
            final AlignData alignment, final byte orientation, final int overlap, final RefGenomeInterface refGenome)
    {
        if(orientation == POS_ORIENT)
        {
            return refGenome.getBaseString(
                    alignment.RefLocation.Chromosome, alignment.RefLocation.end() - overlap + 1, alignment.RefLocation.end());
        }
        else
        {
            return refGenome.getBaseString(
                    alignment.RefLocation.Chromosome, alignment.RefLocation.start(), alignment.RefLocation.start() + overlap - 1);
        }
    }
}
