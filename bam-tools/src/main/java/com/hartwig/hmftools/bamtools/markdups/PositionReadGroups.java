package com.hartwig.hmftools.bamtools.markdups;

import static java.lang.String.format;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.bamtools.ReadGroup;

import htsjdk.samtools.SAMRecord;

public class PositionReadGroups
{
    public final int Position;
    public final Map<String,ReadGroup> ReadGroups;

    public PositionReadGroups(final SAMRecord read)
    {
        ReadGroups = Maps.newHashMap();
        ReadGroups.put(read.getReadName(), new ReadGroup(read));
        Position = read.getAlignmentStart();
    }

    public String toString()
    {
        return format("%d: reads(%d)", Position, ReadGroups.size());
    }
}
