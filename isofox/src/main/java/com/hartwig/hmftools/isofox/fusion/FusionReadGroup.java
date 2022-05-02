package com.hartwig.hmftools.isofox.fusion;

import java.util.List;
import java.util.Map;
import com.google.common.collect.Lists;

public class FusionReadGroup
{
    public final String ReadId;
    public final List<FusionRead> Reads;

    public FusionReadGroup(final String readId, final List<FusionRead> reads)
    {
        ReadId = readId;
        Reads = Lists.newArrayListWithCapacity(reads.size());
        Reads.addAll(reads);
    }

    public int size() { return Reads.size(); }

    public boolean isComplete() { return Reads.size() == 3 || (Reads.size() == 2 && !hasSuppAlignment(Reads)); }

    public boolean hasSuppAlignment() { return hasSuppAlignment(Reads); }

    public static boolean hasSuppAlignment(final List<FusionRead> reads)
    {
        return reads.stream().anyMatch(x -> x.HasSuppAlignment);
    }

    public boolean hasDuplicateRead() { return Reads.stream().anyMatch(x -> x.IsDuplicate); }

    public void merge(final FusionReadGroup other)
    {
        Reads.addAll(other.Reads);
    }

    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", ReadId, Reads.size(), isComplete());
    }

    public String findOtherChromosome(final String chromosome)
    {
        for(FusionRead read : Reads)
        {
            if(!read.MateChromosome.equals(chromosome))
                return read.MateChromosome;

            if(read.SuppData != null)
                return read.SuppData.Chromosome;
        }

        return null;
    }

    public static void mergeChimericReadMaps(
            final Map<String, FusionReadGroup> partialGroups, final List<FusionReadGroup> completeGroups, final Map<String, FusionReadGroup> sourceMap)
    {
        // 1. copies complete groups from the source map into complete groups map
        // 2. checks for a partial match by combining partials and source, and if found removes from partials
        // 3. new partial groups from the source map are copied into the partials map
        // note: source map is logically const
        for(Map.Entry<String, FusionReadGroup> entry : sourceMap.entrySet())
        {
            final FusionReadGroup srcReadGroup = entry.getValue();

            if(srcReadGroup.isComplete())
            {
                completeGroups.add(srcReadGroup);
            }
            else
            {
                // look for an existing incomplete group to add these reads to
                final String readId = entry.getKey();
                FusionReadGroup existingReadGroup = partialGroups.get(readId);

                if(existingReadGroup == null)
                {
                    partialGroups.put(readId, srcReadGroup);
                }
                else
                {
                    existingReadGroup.merge(srcReadGroup);

                    if(existingReadGroup.isComplete())
                    {
                        partialGroups.remove(readId);
                        completeGroups.add(existingReadGroup);
                    }
                }
            }
        }
    }

}
