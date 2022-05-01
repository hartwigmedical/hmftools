package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.isofox.fusion.FusionUtils.suppAlignmentChromosome;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;

public class ReadGroup
{
    public final List<ReadRecord> Reads;

    public ReadGroup(final List<ReadRecord> reads)
    {
        Reads = Lists.newArrayListWithCapacity(reads.size());
        Reads.addAll(reads);
    }

    public ReadGroup(final ReadRecord read)
    {
        Reads = Lists.newArrayListWithCapacity(3);
        Reads.add(read);
    }

    public ReadGroup(final ReadRecord read1, final ReadRecord read2)
    {
        Reads = Lists.newArrayListWithCapacity(3);
        Reads.add(read1);
        Reads.add(read2);
    }

    public final String id() { return Reads.get(0).Id; }

    public int size() { return Reads.size(); }

    public boolean isComplete() { return Reads.size() == 3 || (Reads.size() == 2 && !hasSuppAlignment(Reads)); }

    public boolean hasSuppAlignment() { return hasSuppAlignment(Reads); }

    public static boolean hasSuppAlignment(final List<ReadRecord> reads)
    {
        return reads.stream().anyMatch(x -> x.hasSuppAlignment());
    }

    public boolean hasDuplicateRead() { return Reads.stream().anyMatch(x -> x.isDuplicate()); }

    public void merge(final ReadGroup other)
    {
        Reads.addAll(other.Reads);
    }

    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", id(), Reads.size(), isComplete());
    }

    public String findOtherChromosome(final String chromosome)
    {
        for(ReadRecord read : Reads)
        {
            if(!read.mateChromosome().equals(chromosome))
                return read.mateChromosome();

            if(read.hasSuppAlignment())
                return suppAlignmentChromosome(read.getSuppAlignment());
        }

        return null;
    }

    public static void mergeChimericReadMaps(
            final Map<String,ReadGroup> partialGroups, final List<ReadGroup> completeGroups, final Map<String,ReadGroup> sourceMap)
    {
        // 1. copies complete groups from the source map into complete groups map
        // 2. checks for a partial match by combining partials and source, and if found removes from partials
        // 3. new partial groups from the source map are copied into the partials map
        // note: source map is logically const
        for(Map.Entry<String,ReadGroup> entry : sourceMap.entrySet())
        {
            final ReadGroup srcReadGroup = entry.getValue();

            if(srcReadGroup.isComplete())
            {
                completeGroups.add(srcReadGroup);
            }
            else
            {
                // look for an existing incomplete group to add these reads to
                final String readId = entry.getKey();
                ReadGroup existingReadGroup = partialGroups.get(readId);

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
