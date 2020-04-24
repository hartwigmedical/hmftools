package com.hartwig.hmftools.isofox.common;

public class TransExonRef
{
    public final String GeneId;
    public final int TransId;
    public final String TransName;
    public final int ExonRank;

    public TransExonRef(final String geneId, final int transId, final String transName, final int exonRank)
    {
        GeneId = geneId;
        TransId = transId;
        TransName = transName;
        ExonRank = exonRank;
    }

    public boolean matches(final TransExonRef other)
    {
        return other.TransId == TransId && other.ExonRank == ExonRank;
    }

    public boolean matchesNext(final TransExonRef other)
    {
        // other is one exon ahead of this
        return other.TransId == TransId && other.ExonRank == ExonRank + 1;
    }

    public String toString()
    {
        return String.format("%d:%s:%d", TransId, TransName, ExonRank);
    }

}
