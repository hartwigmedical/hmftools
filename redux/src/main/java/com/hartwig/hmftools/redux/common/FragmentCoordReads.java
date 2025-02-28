package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

public class FragmentCoordReads
{
    public final List<DuplicateGroup> DuplicateGroups;
    public final List<ReadInfo> SingleReads;

    public FragmentCoordReads(final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleReads)
    {
        DuplicateGroups = duplicateGroups != null ? duplicateGroups : Collections.emptyList();
        SingleReads = singleReads != null ? singleReads : Lists.newArrayList(); // may be expanded
    }

    public int duplicateGroupReadCount()
    {
        return DuplicateGroups.stream().mapToInt(x -> x.readCount()).sum();
    }
    public int totalReadCount() { return duplicateGroupReadCount() + SingleReads.size(); }
    public int coordinateCount() { return DuplicateGroups.size() + SingleReads.size(); }

    public int minReadPositionStart()
    {
        int minPositionStart = SingleReads.stream().mapToInt(x -> x.read().getAlignmentStart()).min().orElse(-1);

        for(DuplicateGroup duplicateGroup : DuplicateGroups)
        {
            int minGroupReadPosition = duplicateGroup.reads().stream().mapToInt(x -> x.getAlignmentStart()).min().orElse(-1);

            if(minPositionStart < 0 || minGroupReadPosition > 0 && minGroupReadPosition < minPositionStart)
                minPositionStart = minGroupReadPosition;
        }

        return minPositionStart;
    }

    public String toString()
    {
        return format("duplicateGroups(%d) singles(%d) totalReads(%d)",     DuplicateGroups.size(), SingleReads.size(), totalReadCount());
    }
}
