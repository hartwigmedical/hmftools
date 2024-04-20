package com.hartwig.hmftools.esvee.depth;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;

import htsjdk.samtools.SAMRecord;

class ReadGroup
{
    public final List<SAMRecord> Reads;
    public boolean WaitForAll;

    public ReadGroup(final SAMRecord read, boolean waitForAll)
    {
        Reads = Lists.newArrayList(read);
        WaitForAll = waitForAll;
    }

    public String id() { return Reads.get(0).getReadName(); }

    public String toString()
    {
        return format("id(%s) reads(%d) wait(%s)", id(), Reads.size(), WaitForAll);
    }
}
