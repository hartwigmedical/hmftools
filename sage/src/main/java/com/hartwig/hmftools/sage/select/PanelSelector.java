package com.hartwig.hmftools.sage.select;

import java.util.List;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class PanelSelector<R extends GenomeRegion>
{

    private final List<R> regions;
    private int index = 0;

    public PanelSelector(final List<R> regions)
    {
        this.regions = regions;
    }

    public boolean inPanel(long start, long end)
    {
        if(regions.isEmpty())
        {
            return false;
        }

        R current = current();
        while(index > 0 && current.start() > end)
        {
            index--;
            current = current();
        }

        while(index < regions.size() - 1 && current.end() < start)
        {
            index++;
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
        return regions.get(index);
    }

}
