package com.hartwig.hmftools.telo;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

public class ReadGroup
{
    public final List<SAMRecord> Reads;

    public ReadGroup(final List<SAMRecord> reads)
    {
        Reads = Lists.newArrayListWithCapacity(3);
        Reads.addAll(reads);
    }

    public ReadGroup(final SAMRecord read)
    {
        Reads = Lists.newArrayListWithCapacity(3);
        Reads.add(read);
    }

    public final String id() { return Reads.get(0).getReadName(); }

    public boolean isComplete() { return Reads.size() == 2 || !Reads.get(0).getReadPairedFlag(); }

    public void merge(final ReadGroup other)
    {
        Reads.addAll(other.Reads);
    }

    public String toString() { return String.format("%s reads(%d) complete(%s)", id(), Reads.size(), isComplete()); }

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
