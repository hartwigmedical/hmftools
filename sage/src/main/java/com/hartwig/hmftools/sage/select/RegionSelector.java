package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.Optional;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class RegionSelector<R extends GenomeRegion>
{
    private final List<R> mRegions;
    private int mIndex = 0;

    public RegionSelector(final List<R> regions)
    {
        this.mRegions = regions;
    }

    @NotNull
    public Optional<R> select(long position)
    {
        if(mRegions.isEmpty())
        {
            return Optional.empty();
        }

        R current = current();
        while(mIndex > 0 && current.start() > position)
        {
            mIndex--;
            current = current();
        }

        while(mIndex < mRegions.size() - 1 && current.end() < position)
        {
            mIndex++;
            current = current();
        }

        if(position >= current.start() && position <= current.end())
        {
            return Optional.of(current);
        }

        return Optional.empty();
    }

    @NotNull
    private R current()
    {
        return mRegions.get(mIndex);
    }

}
