package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.Optional;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class RegionSelector<R extends GenomeRegion>
{

    private final List<R> regions;
    private int index = 0;

    public RegionSelector(final List<R> regions)
    {
        this.regions = regions;
    }

    @NotNull
    public Optional<R> select(long position)
    {
        if(regions.isEmpty())
        {
            return Optional.empty();
        }

        R current = current();
        while(index > 0 && current.start() > position)
        {
            index--;
            current = current();
        }

        while(index < regions.size() - 1 && current.end() < position)
        {
            index++;
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
        return regions.get(index);
    }

}
