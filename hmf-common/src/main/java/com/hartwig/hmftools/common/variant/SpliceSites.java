package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.genome.region.Strand.FORWARD;

import com.hartwig.hmftools.common.genome.region.Strand;

public final class SpliceSites
{
    public static int getDonorPosition(int position, int exonEnd, Strand strand)
    {
        // D-1 refs to the last base of the exon, D+1 the first base after, ie first of the intron, and there is no D-zero
        if(strand == FORWARD)
        {
            return position > exonEnd ? position - exonEnd : position - exonEnd - 1;
        }
        else
        {
            return position < exonEnd ? exonEnd - position : exonEnd - position - 1;
        }
    }

    public static int getAcceptorPosition(int position, int exonEnd, Strand strand)
    {
        // A-1 refs to the first base of the exon, A+1 the first base before, ie first preceding in the intron, and there is no A-zero
        if(strand == FORWARD)
        {
            return position < exonEnd ? exonEnd - position : exonEnd - position - 1;
        }
        else
        {
            return position > exonEnd ? position - exonEnd : position - exonEnd - 1;
        }
    }
}
