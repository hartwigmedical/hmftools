package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.ReadGroup;

import htsjdk.samtools.SAMRecord;

public class PositionReadGroups
{
    public final int Position;
    public final List<ReadGroup> ReadGroups;

    public PositionReadGroups(final SAMRecord record)
    {
        ReadGroups = Lists.newArrayList();
        ReadGroups.add(new ReadGroup(record));
        Position = record.getAlignmentStart();
    }

    public String toString()
    {
        return format("%d: reads(%d)", Position, ReadGroups.size());
    }
}
