package com.hartwig.hmftools.isofox.fusion;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;

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

    public void merge(final ReadGroup other)
    {
        Reads.addAll(other.Reads);
    }

    public static void mergeChimericReadMaps(final Map<String,ReadGroup> destMap, final Map<String,ReadGroup> sourceMap)
    {
        for(Map.Entry<String,ReadGroup> entry : sourceMap.entrySet())
        {
            ReadGroup readsById = destMap.get(entry.getKey());

            if(readsById == null)
            {
                destMap.put(entry.getKey(), entry.getValue());
            }
            else
            {
                readsById.merge(entry.getValue());
            }
        }
    }

}
