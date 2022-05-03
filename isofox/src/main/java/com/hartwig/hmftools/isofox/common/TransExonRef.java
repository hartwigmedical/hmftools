package com.hartwig.hmftools.isofox.common;

import java.util.List;

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

    private boolean matches(final TransExonRef other, int maxDiffVsOther)
    {
        // if argument is +ve then the other trans ref's exon must be equal or higher than this one by the permitted diff
        if(TransId != other.TransId)
            return false;

        if(maxDiffVsOther >= 0)
            return other.ExonRank >= ExonRank && other.ExonRank <= ExonRank + maxDiffVsOther;
        else
            return other.ExonRank >= ExonRank + maxDiffVsOther && other.ExonRank <= ExonRank;
    }

    public static boolean hasTranscriptExonMatch(final List<TransExonRef> list1, final List<TransExonRef> list2)
    {
        return list1.stream().anyMatch(x -> list2.stream().anyMatch(y -> x.matches(y)));
    }

    public static void mergeUnique(final List<TransExonRef> existingRefs, final List<TransExonRef> newRefs)
    {
        for(TransExonRef transExonRef : newRefs)
        {
            if(existingRefs.stream().noneMatch(x -> x.TransId == transExonRef.TransId && x.ExonRank == transExonRef.ExonRank))
            {
                existingRefs.add(transExonRef);
            }
        }

    }
    public static boolean hasTranscriptExonMatch(final List<TransExonRef> list1, final List<TransExonRef> list2, int maxList2Diff)
    {
        // list 2 is allowed a buffer of difference
        return list1.stream().anyMatch(x -> list2.stream().anyMatch(y -> x.matches(y, maxList2Diff)));
    }



    public String toString()
    {
        return String.format("%d:%s:%d", TransId, TransName, ExonRank);
    }

}
