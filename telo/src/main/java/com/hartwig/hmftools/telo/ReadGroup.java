package com.hartwig.hmftools.telo;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class ReadGroup
{
    public final List<ReadRecord> Reads;

    public ReadGroup(final List<ReadRecord> reads)
    {
        Reads = Lists.newArrayListWithCapacity(3);
        Reads.addAll(reads);
    }

    public ReadGroup(final ReadRecord read)
    {
        Reads = Lists.newArrayListWithCapacity(3);
        Reads.add(read);
    }

    public final String id() { return Reads.get(0).Id; }

    public boolean isComplete() { return Reads.size() == 3 || (Reads.size() == 2 && !hasSuppAlignment(Reads)); }

    public static boolean hasSuppAlignment(final List<ReadRecord> reads)
    {
        return reads.stream().anyMatch(x -> x.hasSuppAlignment());
    }

    public void merge(final ReadGroup other)
    {
        Reads.addAll(other.Reads);
    }

    public String toString() { return String.format("%s reads(%d) complete(%s)", id(), Reads.size(), isComplete()); }

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

    public static final String SUPP_ALIGNMENT_DELIM = ",";

    public static String suppAlignmentChromosome(final String suppAlignment)
    {
        if(suppAlignment == null)
            return null;

        final String[] items = suppAlignment.split(SUPP_ALIGNMENT_DELIM);
        return items.length >= 5 ? items[0] : null;
    }

    public static void mergeReadGroups(
            final Map<String,ReadGroup> partialGroups, final List<ReadGroup> completeGroups, final Map<String,ReadGroup> sourceMap)
    {
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
