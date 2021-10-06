package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.pave.PaveConstants.GENE_UPSTREAM_DISTANCE;

import com.hartwig.hmftools.common.gene.TranscriptData;

public final class PaveUtils
{
    public static boolean withinTransRange(final TranscriptData transData, int posStart, int posEnd)
    {
        if(transData.posStrand())
            return positionsOverlap(posStart, posEnd, transData.TransStart - GENE_UPSTREAM_DISTANCE, transData.TransEnd);
        else
            return positionsOverlap(posStart, posEnd, transData.TransStart, transData.TransEnd + GENE_UPSTREAM_DISTANCE);
    }

}
