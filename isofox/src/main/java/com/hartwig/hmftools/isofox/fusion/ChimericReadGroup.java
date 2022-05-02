package com.hartwig.hmftools.isofox.fusion;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.isofox.common.ReadRecord;

public class ChimericReadGroup
{
    public final List<ReadRecord> Reads;

    public ChimericReadGroup(final ReadRecord read)
    {
        Reads = Lists.newArrayListWithCapacity(2);
        Reads.add(read);
    }

    public ChimericReadGroup(final ReadRecord read1, final ReadRecord read2)
    {
        Reads = Lists.newArrayListWithCapacity(2);
        Reads.add(read1);
        Reads.add(read2);
    }

    public final String id() { return Reads.get(0).Id; }

    public int size() { return Reads.size(); }

    public boolean isComplete() { return Reads.size() == 3 || (Reads.size() == 2 && !hasSuppAlignment()); }

    public boolean hasSuppAlignment() { return Reads.stream().anyMatch(x -> x.hasSuppAlignment()); }

    public String toString()
    {
        return String.format("%s reads(%d) complete(%s)", id(), Reads.size(), isComplete());
    }
}
