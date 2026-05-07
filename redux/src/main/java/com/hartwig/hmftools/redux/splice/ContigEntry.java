package com.hartwig.hmftools.redux.splice;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;

public class ContigEntry
{
    public final String ContigName;
    public final String GeneId;
    public final String GeneName;
    public final String TransName;
    public final String Chromosome;

    // exon spans on the forward genomic strand, sorted ascending. Each span is 1-based, inclusive.
    public final List<BaseRegion> ExonSpans;

    public ContigEntry(
            final String contigName, final String geneId, final String geneName, final String transName,
            final String chromosome, final List<BaseRegion> exonSpans)
    {
        ContigName = contigName;
        GeneId = geneId;
        GeneName = geneName;
        TransName = transName;
        Chromosome = chromosome;
        ExonSpans = exonSpans;
    }

    public int contigLength()
    {
        int length = 0;
        for(BaseRegion span : ExonSpans)
            length += span.baseLength();
        return length;
    }
}
