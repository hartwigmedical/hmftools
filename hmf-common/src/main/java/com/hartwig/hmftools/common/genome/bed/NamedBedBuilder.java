package com.hartwig.hmftools.common.genome.bed;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class NamedBedBuilder implements Consumer<NamedBed>
{
    private final Map<String, NamedBedBuilderChromosome> chromosomeMap;

    public NamedBedBuilder()
    {
        this.chromosomeMap = Maps.newHashMap();
    }

    public boolean addBed(NamedBed bed)
    {
        return chromosomeMap.computeIfAbsent(bed.chromosome(), NamedBedBuilderChromosome::new).addBed(bed);
    }

    @NotNull
    public List<NamedBed> build()
    {
        List<NamedBed> result = Lists.newArrayList();
        for(String contig : chromosomeMap.keySet())
        {
            result.addAll(chromosomeMap.get(contig).build());
        }

        Collections.sort(result);
        return result;
    }

    @Override
    public void accept(final NamedBed namedBed)
    {
        addBed(namedBed);
    }

    static final class NamedBedBuilderChromosome
    {
        private final String chromosome;
        @NotNull
        private final List<NamedBed> regions;

        public NamedBedBuilderChromosome(@NotNull final String chromosome)
        {
            this.chromosome = chromosome;
            this.regions = Lists.newArrayList();
        }

        @NotNull
        public List<NamedBed> build()
        {
            return regions;
        }

        public boolean addBed(@NotNull final NamedBed bed)
        {
            if(regions.isEmpty())
            {
                regions.add(bed);
                return true;
            }

            int regionBeforeStart = indexOfRegionEndingBefore(bed.start());

            // Add to end
            if(regionBeforeStart == regions.size() - 1)
            {
                regions.add(bed);
                return true;
            }

            int regionAfterEnd = indexOfRegionStartingAfter(bed.end());

            // Add to start
            if(regionAfterEnd == 0)
            {
                regions.add(0, bed);
                return true;
            }

            // Add in between
            if(regionAfterEnd == regionBeforeStart + 1)
            {
                regions.add(regionAfterEnd, bed);
                return true;
            }

            int indexOfStart = indexOf(bed.start());
            int indexOfEnd = indexOf(bed.end());
            if(indexOfStart == indexOfEnd)
            {
                return false;
            }

            throw new IllegalArgumentException("Don't currently handle partial overlaps");
        }

        private int indexOfRegionStartingAfter(long position)
        {
            for(int i = 0; i < regions.size(); i++)
            {
                GenomeRegion region = regions.get(i);
                if(region.start() > position)
                {
                    return i;
                }
            }

            return -1;
        }

        private int indexOfRegionEndingBefore(long position)
        {
            int result = -1;
            for(int i = regions.size() - 1; i >= 0; i--)
            {
                GenomeRegion region = regions.get(i);
                if(region.end() < position)
                {
                    return i;
                }
            }

            return result;
        }

        private int indexOf(long position)
        {
            for(int i = regions.size() - 1; i >= 0; i--)
            {
                GenomeRegion region = regions.get(i);
                if(position >= region.start() && position <= region.end())
                {
                    return i;
                }
            }

            return -1;
        }
    }
}
