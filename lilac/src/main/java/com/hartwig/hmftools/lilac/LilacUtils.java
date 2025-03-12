package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import java.util.Arrays;
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

    public static List<Integer> arrayToList(final int[] array)
    {
        List<Integer> list = Lists.newArrayListWithExpectedSize(array.length);
        Arrays.stream(array).forEach(x -> list.add(x));
        return list;
    }

    public static List<String> arrayToList(final char[] array)
    {
        List<String> list = Lists.newArrayListWithExpectedSize(array.length);

        for(int i = 0; i < array.length; ++i)
            list.add(String.valueOf(array[i]));

        return list;
    }

    public static int listMin(final List<Integer> list)
    {
        return list.stream().mapToInt(x -> x).min().orElse(0);
    }

    public static int listMax(final List<Integer> list)
    {
        return list.stream().mapToInt(x -> x).max().orElse(0);
    }

    public static boolean namesMatch(final Set<String> list1, final Set<String> list2)
    {
        if(list1.size() != list2.size())
            return false;

        return !list1.stream().noneMatch(x -> list2.contains(x));
    }

}

