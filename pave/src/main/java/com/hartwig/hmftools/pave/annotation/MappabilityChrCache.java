package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.common.variant.PaveVcfTags.MAPPABILITY;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.pave.VariantData;

public class MappabilityChrCache
{
    public final String Chromosome;

    private int mCurrentIndex;
    private final List<MapEntry> mEntries;
    private boolean mComplete;

    public MappabilityChrCache(final String chromosome)
    {
        Chromosome = chromosome;
        mCurrentIndex = 0;
        mEntries = Lists.newArrayList();
        mComplete = false;
    }

    public void addEntry(final int posStart, final int posEnd, final double mappability)
    {
        mEntries.add(new MapEntry(new BaseRegion(posStart, posEnd), mappability));
    }

    public boolean isComplete() { return mComplete; }
    public void setComplete() { mComplete = true; }
    public void clear() { mEntries.clear(); }
    public int entryCount() { return mEntries.size(); }

    public void annotateVariant(final VariantData variant)
    {
        if(mEntries.isEmpty() || mEntries.get(mEntries.size() - 1).Region.end() < variant.Position)
            return;

        int position = variant.Position;

        for(; mCurrentIndex < mEntries.size(); ++mCurrentIndex)
        {
            MapEntry entry = mEntries.get(mCurrentIndex);

            if(entry.Region.end() < position)
                continue;

            if(entry.Region.containsPosition(variant.Position))
            {
                setMappability(variant, entry.Mappability);
            }
            else if(variant.Position < entry.Region.start() && mCurrentIndex > 0)
            {
                // take previous if the next is past this variant
                --mCurrentIndex;
                MapEntry prevEntry = mEntries.get(mCurrentIndex);
                setMappability(variant, prevEntry.Mappability);
            }

            return;
        }
    }

    private void setMappability(final VariantData variant, double mappability)
    {
        if(!variant.context().getCommonInfo().hasAttribute(MAPPABILITY))
        {
            variant.context().getCommonInfo().putAttribute(MAPPABILITY, mappability);
        }
    }

    private class MapEntry
    {
        public final BaseRegion Region;
        public final double Mappability;

        public MapEntry(final BaseRegion region, final double mappability)
        {
            Region = region;
            Mappability = mappability;
        }

        public String toString() { return String.format("%s map(%.4f)", Region, Mappability); }
    }
}
