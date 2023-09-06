package com.hartwig.hmftools.pave.annotation;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.pave.VariantData;

public class PonChrCache
{
    public final String Chromosome;

    private final Map<Integer,List<PonVariantData>> mPositionMap;
    private boolean mComplete;

    public PonChrCache(final String chromosome)
    {
        Chromosome = chromosome;
        mPositionMap = Maps.newHashMap();
        mComplete = false;
    }

    public void addEntry(
            final int position, final String ref, final String alt, final int samples, final int maxSampleReads, final int totalSampleReads)
    {
        List<PonVariantData> posList = mPositionMap.get(position);

        if(posList == null)
        {
            posList = Lists.newArrayList();
            mPositionMap.put(position, posList);
        }

        posList.add(new PonVariantData(ref, alt, samples, maxSampleReads, totalSampleReads));
    }

    public boolean isComplete() { return mComplete; }
    public void setComplete() { mComplete = true; }
    public void clear() { mPositionMap.clear(); }
    public int entryCount() { return mPositionMap.values().stream().mapToInt(x -> x.size()).sum(); }

    public String toString() { return format("chr(%s) entries(%d)", Chromosome, mPositionMap.size()); }

    public PonVariantData getPonData(final VariantData variant)
    {
        return getPonData(variant.Position, variant.Ref, variant.Alt);
    }

    public boolean hasEntry(final int position, final String ref, final String alt)
    {
        return getPonData(position, ref, alt) != null;
    }

    private PonVariantData getPonData(final int position, final String ref, final String alt)
    {
        List<PonVariantData> posList = mPositionMap.get(position);

        if(posList == null)
            return null;

        return posList.stream().filter(x -> x.matches(ref, alt)).findFirst().orElse(null);
    }
}
