package com.hartwig.hmftools.geneutils.mapping;

import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;

public class LiftOverRegion
{
    public final String GeneId;
    public final int GeneStart;
    public final int GeneEnd;

    public LiftOverRegion(final String geneId, final int geneStart, final int geneEnd)
    {
        GeneId = geneId;
        GeneStart = geneStart;
        GeneEnd = geneEnd;
    }

    public boolean positionMatches(int start, int end, int buffer)
    {
        if(start < GeneStart - buffer || start > GeneStart + buffer)
            return false;

        if(end < GeneEnd - buffer || end > GeneEnd + buffer)
            return false;

        return true;
    }
}
