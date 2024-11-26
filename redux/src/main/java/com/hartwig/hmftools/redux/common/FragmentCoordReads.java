package com.hartwig.hmftools.redux.common;

import static java.lang.String.format;

import java.util.Collections;
import java.util.List;

public class FragmentCoordReads
{
    public final List<DuplicateGroup> DuplicateGroups;
    public final List<ReadInfo> SingleReads;

    public FragmentCoordReads(final List<DuplicateGroup> duplicateGroups, final List<ReadInfo> singleReads)
    {
        DuplicateGroups = duplicateGroups != null ? duplicateGroups : Collections.emptyList();
        SingleReads = singleReads != null ? singleReads : Collections.emptyList();
    }

    public int totalReadCount()
    {
        return DuplicateGroups.stream().mapToInt(x -> x.readCount()).sum() + SingleReads.size();
    }

    public int coordinateCount()
    {
        return DuplicateGroups.size() + SingleReads.size();
    }

    public String toString()
    {
        return format("duplicateGroups(%d) singles(%d) totalReads(%d)", DuplicateGroups.size(), SingleReads.size(), totalReadCount());
    }
}
