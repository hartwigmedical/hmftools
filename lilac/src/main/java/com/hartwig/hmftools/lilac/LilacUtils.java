package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptData;

public class LilacUtils
{
    public static int calcNucelotideLocus(final List<TranscriptData> transcripts, int position)
    {
        for(TranscriptData transData : transcripts)
        {
            if(!positionWithin(position, transData.CodingStart, transData.CodingEnd))
                continue;

            // locus is a zero-based index, so the first coding base has locus of 0
            return calcCodingBases(transData, position).CodingBases - 1;
        }
        return -1;
    }

    public static List<Integer> formRange(int start, int end)
    {
        List<Integer> indices = Lists.newArrayList();
        for(int i = start; i <= end; ++i)
        {
            indices.add(i);
        }

        return indices;

    }

    public static int listMin(final List<Integer> list)
    {
        return list.stream().mapToInt(x -> x).min().orElse(0);
    }

    public static int listMax(final List<Integer> list)
    {
        return list.stream().mapToInt(x -> x).max().orElse(0);
    }

    public static <T> boolean namesMatch(final Set<T> list1, final Set<T> list2)
    {
        if(list1.size() != list2.size())
            return false;

        return !list1.stream().noneMatch(x -> list2.contains(x));
    }

}

