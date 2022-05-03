package com.hartwig.hmftools.isofox.fusion;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.isofox.common.TransExonRef;

public class FusionTransExon
{
    public final int TransId;
    public final int ExonRank;

    public FusionTransExon(final int transId, final int exonRank)
    {
        TransId = transId;
        ExonRank = exonRank;
    }

    public FusionTransExon(final TransExonRef transExonRef)
    {
        TransId = transExonRef.TransId;
        ExonRank = transExonRef.ExonRank;
    }

    public static List<FusionTransExon> fromList(final List<TransExonRef> transExonRefs)
    {
        return transExonRefs.stream().map(x -> new FusionTransExon(x)).collect(Collectors.toList());
    }

    public boolean matches(final FusionTransExon other)
    {
        return other.TransId == TransId && other.ExonRank == ExonRank;
    }

    public static boolean matches(final TransExonRef teRef, final FusionTransExon fusionTeRef)
    {
        return teRef.TransId == fusionTeRef.TransId && teRef.ExonRank == fusionTeRef.ExonRank;
    }

    public static boolean hasMatchWithinRange(final List<TransExonRef> teRefs, final List<FusionTransExon> fusionTeRefs, int maxList2Diff)
    {
        // list 2 is allowed a buffer of difference
        return teRefs.stream().anyMatch(x -> fusionTeRefs.stream().anyMatch(y -> matchesWithinRange(x, y, maxList2Diff)));
    }

    private static boolean matchesWithinRange(final TransExonRef teRef, final FusionTransExon fusionTeRef, int maxDiffVsOther)
    {
        // if argument is +ve then the other trans ref's exon must be equal or higher than this one by the permitted diff
        if(teRef.TransId != fusionTeRef.TransId)
            return false;

        if(maxDiffVsOther >= 0)
            return fusionTeRef.ExonRank >= teRef.ExonRank && fusionTeRef.ExonRank <= teRef.ExonRank + maxDiffVsOther;
        else
            return fusionTeRef.ExonRank >= teRef.ExonRank + maxDiffVsOther && fusionTeRef.ExonRank <= teRef.ExonRank;
    }

    public static void mergeUnique(final List<FusionTransExon> existingRefs, final List<FusionTransExon> newRefs)
    {
        for(FusionTransExon transExonRef : newRefs)
        {
            if(existingRefs.stream().noneMatch(x -> x.matches(transExonRef)))
            {
                existingRefs.add(transExonRef);
            }
        }
    }

    public String toString()
    {
        return String.format("%d:%d", TransId, ExonRank);
    }
}
