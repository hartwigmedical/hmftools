package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class PanelSelector<R extends GenomeRegion>
{
    private final List<R> mRegions;
    private int mIndex = 0;

    public PanelSelector(final List<R> regions)
    {
        this.mRegions = regions;
    }

    public boolean inPanel(long start, long end)
    {
        if(mRegions.isEmpty())
        {
            return false;
        }

        R current = current();
        while(mIndex > 0 && current.start() > end)
        {
            mIndex--;
            current = current();
        }

        while(mIndex < mRegions.size() - 1 && current.end() < start)
        {
            mIndex++;
            current = current();
        }

        if(start <= current.end() && end >= current.start())
        {
            return true;
        }

        return false;
    }

    @NotNull
    private R current()
    {
        return mRegions.get(mIndex);
    }

}
